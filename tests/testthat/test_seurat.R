context("Creating .h5mu files from Seurat objects")
library(Seurat)
library(MuDataSeurat)
library(hdf5r)
library(fs)  # for file_temp()

fileh5mu <- paste0(file_temp(), ".h5mu")
fileh5ad <- paste0(file_temp(), ".h5ad")

# Hard-coded value to inject
true_val <- 0.1234569
true_val_i <- 3
true_val_j <- 7

test_that("a Seurat object with two modalities can be created, written, and read back", {
    x <- matrix(rnorm(1000), ncol = 100)
    y <- matrix(rnorm(2000), ncol = 100)
    obs_names <- paste("obs", 1:100, sep = "-")

    colnames(x) <- obs_names
    colnames(y) <- obs_names

    rownames(x) <- paste("x-var", 1:10, sep = "-")
    rownames(y) <- paste("y-var", 1:20, sep = "-")

    # Inject a value to be checked later
    x[true_val_i, true_val_j] <- true_val

    assay_x <- CreateAssayObject(counts = x)
    assay_y <- CreateAssayObject(counts = y)

    srt <- CreateSeuratObject(assay_x, assay = "x")
    srt[["y"]] <- assay_y
    
    # Writing
    outfile <- fileh5mu
    result <- WriteH5MU(srt, outfile)

    # Assert the data is saved
    expect_true(result)

    # Read back
    h5 <- H5File$new(outfile, mode = "r")

    # Check all the assays are written
    assays <- h5[['mod']]$names
    assays_orig <- sort(names(srt))
    expect_equal(assays, assays_orig)

    h5$close_all()
})

test_that("a Seurat object can be created from an .h5mu file", {
    srt <- ReadH5MU(fileh5mu)
    expect_equal(names(srt)[1], "x")
    expect_equal(names(srt)[2], "y")
})

test_that("values in matrices are properly recovered", {
    srt <- ReadH5MU(fileh5mu)
    expect_equal(srt[["x"]]@counts[true_val_i, true_val_j], true_val)
    expect_equal(srt[["x"]]@data[true_val_i, true_val_j], true_val)
})

test_that(".h5ad I/O works", {
    x <- matrix(rnorm(200), ncol = 20)
    colnames(x) <- paste("obs", 1:20, sep = "-")
    rownames(x) <- paste("var", 1:10, sep = "-")
    x[true_val_i, true_val_j] <- true_val
    srt <- CreateSeuratObject(counts = x)

    # Write
    result <- WriteH5AD(srt, fileh5ad)
    expect_true(result)

    # Read back
    srt <- ReadH5AD(fileh5ad)
    expect_equal(srt@assays$RNA@counts[true_val_i, true_val_j], true_val)
})