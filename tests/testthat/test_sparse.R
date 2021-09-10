context("Creating .h5ad and .h5mu files with sparse matrices")
library(Seurat)
library(MuDataSeurat)
library(Matrix)
library(hdf5r)
library(fs)  # for file_temp()

# csc
fileh5mu_r <- paste0(file_temp(), ".h5mu")
fileh5ad_r <- paste0(file_temp(), ".h5ad")
# csc
fileh5mu_c <- paste0(file_temp(), ".h5mu")
fileh5ad_c <- paste0(file_temp(), ".h5ad")


nobs <- 10
nvar <- 20
nvar2 <- 31

obs_names <- paste("obs", 1:nobs, sep="-")
var_names <- paste("var", 1:nvar, sep="-")

# Hard-coded value to inject
true_val <- 0.1234569
true_val_i <- 3
true_val_j <- 7

test_that("dgCMatrix can be written to .h5ad", {
    x <- rnbinom(n = nobs * nvar, prob = .95, size = 10)
    x <- Matrix(matrix(x, ncol = nobs), sparse = TRUE)  # => dgCMatrix
    x[true_val_i,true_val_j] <- true_val

    colnames(x) <- obs_names
    rownames(x) <- var_names

    srt <- CreateSeuratObject(counts = x)

    expect_true(WriteH5AD(srt, fileh5ad_c))

    h5 <- H5File$new(fileh5ad_c, mode="a")
    expect_equal(h5attr(h5[["X"]], "encoding-type"), "csc_matrix")
    h5$close()
})

test_that("dgRMatrix can be written to .h5ad", {
    # Seurat objects would not accept dgRMatrix,
    # and there is no native dgRMatrix -> dgCMatrix conversion.
    # We will write a transposed matrix as CSC and then change it to CSR manually
    # to emulate CSR matrices written by other tools.

    x <- rnbinom(n = nobs * nvar, prob = .95, size = 10)
    x <- x[x != 0]
    i <- sample(1:nvar, length(x), replace = T)
    j <- sample(1:nobs, length(x), replace = T)
    
    x_c <- sparseMatrix(i = i, j = j, x = x, dims = c(nvar, nobs), repr = "C")    # => dgCMatrix
    x_c[true_val_i,true_val_j] <- true_val

    x_r <- sparseMatrix(i = i, j = j, x = x, dims = c(nvar, nobs), repr = "R")  # => dgRMatrix
    # Value assignment will automatically convert it to dgTMatrix,
    # and there is no dgTMatrix -> dgRMatrix coercion,
    # so in this case we won't inject the value.
    # x_r[true_val_i,true_val_j] <- true_val

    colnames(x_c) <- obs_names
    rownames(x_c) <- var_names

    srt <- CreateSeuratObject(counts = x_c)

    expect_true(WriteH5AD(srt, fileh5ad_r))

    h5 <- H5File$new(fileh5ad_r, mode="a")
    h5x <- h5[["X"]]
    expect_equal(h5attr(h5x, "encoding-type"), "csc_matrix")

    h5x$link_delete("indices")
    h5x[["indices"]] <- x_r@j
    h5x$link_delete("indptr")
    h5x[["indptr"]] <- x_r@p
    h5x$link_delete("data")
    h5x[["data"]] <- x_r@x
    h5attr(h5x, "encoding-type") <- "csr_matrix"

    h5$close()
})

test_that("dgCMatrix can be read from .h5ad", {
    srt <- ReadH5AD(fileh5ad_c)

    counts <- srt@assays[[1]]@counts
    expect_true("dgCMatrix" %in% class(counts))
    expect_equal(dim(counts), c(nvar, nobs))
    expect_equal(counts[true_val_i, true_val_j], true_val)
    expect_equal(rownames(srt), var_names)
    expect_equal(colnames(srt), obs_names)
})

test_that("dgRMatrix can be read from .h5ad", {
    srt <- ReadH5AD(fileh5ad_r)
    counts <- srt@assays[[1]]@counts
    
    # Seurat only support dgCMatrix as counts
    expect_true("dgCMatrix" %in% class(counts))
    expect_equal(dim(counts), c(nvar, nobs))
    
    expect_equal(rownames(srt), var_names)
    expect_equal(colnames(srt), obs_names)
})


