#' @importFrom hdf5r is_hdf5 H5File
open_and_check_mudata <- function(filename) {
    if (readChar(filename, 6) != "MuData") {
        if (is_hdf5(filename)) {
            warning("The HDF5 file was not created by MuData tooling, we can't guarantee that everything will work correctly", call.=FALSE)
        } else (
            stop("The file is not an HDF5 file", call.=FALSE)
        )
    }
    H5File$new(filename, mode="r")
}

#' @import hdf5r
open_anndata <- function(filename) {
  # PATH/filename.h5mu/mod/rna => read a single modality from the .h5mu file
  path_fragments <- strsplit(filename, "\\.h5mu")[[1]]
  if (length(path_fragments) == 1) {
    h5 <- H5File$new(filename, mode="r")
  } else {
    h5 <- H5File$new(paste0(path_fragments[1], ".h5mu"), mode="r")
    mod_path <- path_fragments[2]
    if (substr(mod_path, 1, 4) != "/mod") {
      mod_path <- paste0("/mod", mod_path)
    }
    h5 <- h5[[mod_path]]
  }
  h5
}

missing_on_read <- function(loc, desc = "") {
  details <- ""
  if (!is.null(desc) && desc != "") {
    details <- paste0("Seurat does not support ", desc, ".")
  }
  warning(paste0("Missing on read: ", loc, ". ", details))
}

read_table_encv1 <- function(dataset, set_index = TRUE) {
  columns <- names(dataset)
  columns <- columns[columns != "__categories"]

  col_list <- lapply(columns, function(name) {
    values <- dataset[[name]]$read()
    values_attr <- tryCatch({
      h5attributes(dataset[[name]])
    }, error = function(e) {
      list()
    })
    if (length(values_attr) > 0) {
      if ("categories" %in% names(values_attr)) {
        # Make factors out of categorical data
        ref <- values_attr$categories
        values_labels <- ref$dereference(obj = NULL)[[1]]
        # NOTE: number of labels have to be strictly matching the number of unique integer values.
        values_notna <- unique(values)
        values_notna <- values_notna[!is.na(values_notna)]
        values <- factor(as.integer(values), labels = values_labels$read()[1:length(values_notna)])
      }
    }
    values
  })
  table <- data.frame(Reduce(cbind.data.frame, col_list))
  colnames(table) <- columns
  table
}

read_column <- function(column, etype, eversion) {
  values <- NULL
  if (etype == "categorical") {
    if (eversion == "0.2.0") {
      codes <- column[["codes"]]$read()
      categories <- column[["categories"]]$read()

      # NOTE: number of labels have to be strictly matching the number of unique integer values.
      codes_notna <- unique(codes)
      codes_notna <- codes_notna[!is.na(codes_notna)]

      values <- factor(as.integer(codes), labels = categories[1:length(codes_notna)])
    } else {
      warning(paste0("Cannot recognise encoving-version ", eversion))
    }
  } else {
    values <- column$read()
  }
  # } else {
  #   stop(paste0("Cannot recognise encoving-type ", etype))
  # }
  values
}

read_table_encv2 <- function(dataset, set_index = TRUE) {
  columns <- names(dataset)

  col_list <- lapply(columns, function(name) {

    col_attr <- tryCatch({
      h5attributes(dataset[[name]])
    }, error = function(e) {
      list("encoding-type" = NULL)
    })

    values <- read_column(dataset[[name]], col_attr$`encoding-type`, col_attr$`encoding-version`)

    values
  })
  table <- data.frame(Reduce(cbind.data.frame, col_list))
  colnames(table) <- columns
  table
}

read_table <- function(dataset, set_index = TRUE) {
  if ("H5Group" %in% class(dataset)) {
    # Table is saved as a group rather than a dataset
    dataset_attr <- tryCatch({
      h5attributes(dataset)
    }, error = function(e) {
      list("_index" = "_index")
    })
    indexcol <- "_index"
    if ("_index" %in% names(dataset_attr)) {
      indexcol <- dataset_attr$`_index`
    }

    encv <- "0.1.0"  # some encoding version by default
    if ("encoding-version" %in% names(dataset_attr)) {
      encv <- dataset_attr$`encoding-version`
    }

    if (encv == "0.1.0") {
      table <- read_table_encv1(dataset, set_index)
    } else if (encv == "0.2.0") {
      table <- read_table_encv2(dataset, set_index)
    } else {
      stop(paste0("Encoding version ", encv, " is not recognised."))
    }

    columns <- colnames(table)

    if ((indexcol %in% colnames(table)) && set_index) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }

    # Fix column order
    if ("column-order" %in% names(dataset_attr)) {
      ordered_columns <- dataset_attr[["column-order"]]
      # Do not consider index as a column
      ordered_columns <- ordered_columns[ordered_columns != indexcol]
      table <- table[,ordered_columns[ordered_columns %in% columns],drop=FALSE]
    }
  } else {
    table <- dataset$read()
    dataset_attr <- h5attributes(dataset)

    indexcol <- "_index"
    if ("_index" %in% names(dataset_attr)) {
      indexcol <- dataset_attr$`_index`
    }

    if ((indexcol %in% colnames(table)) && set_index) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }
  }
  table
}

#' @import Matrix
read_matrix <- function(dataset) {
  if ("data" %in% names(dataset) && "indices" %in% names(dataset) && "indptr" %in% names(dataset)) {
      i <- dataset[["indices"]]$read()
      p <- dataset[["indptr"]]$read()
      x <- dataset[["data"]]$read()

      rowwise <- FALSE
      if ("encoding-type" %in% h5attr_names(dataset)) {
        rowwise <- h5attr(dataset, "encoding-type") == "csr_matrix"
      }

      if ("shape" %in% h5attr_names(dataset)) {
        X_dims <- h5attr(dataset, "shape")
      } else {
        X_dims <- c(length(p) - 1, max(i) + 1)
        if (rowwise) {
          X_dims <- rev(X_dims)
        }
      }

      # X is a dgCMatrix.
      # No direct dgCMatrix -> dgRMatrix coersion provided in Matrix.
      if (rowwise) {
        X <- Matrix::Matrix(0, X_dims[2], X_dims[1], doDiag = FALSE)
      } else {
        X <- Matrix::Matrix(0, X_dims[1], X_dims[2], doDiag = FALSE)
      }
      X@i <- i
      X@p <- p
      X@x <- x

      if (rowwise) {
        X
      } else {
        Matrix::t(X)
      }
    } else {
      dataset$read()
    }
}

#' @import Matrix
read_layers_to_assay <- function(root, modalityname="") {
  X <- read_matrix(root[['X']])

  var <- read_table(root[['var']])
  if (any(grepl("_", rownames(var)))) {
    example_which <- grep("_", rownames(var))[1]
    example_before <- rownames(var)[example_which]
    rownames(var) <- gsub("_", "-", rownames(var))
    example_after <- rownames(var)[example_which]
    warning(paste0("The var_names from modality ", modalityname, " have been renamed as feature names cannot contain '_'.",
      " E.g. ", example_before, " -> ", example_after, "."))
  }

  obs <- read_table(root[['obs']])
  if (is("obs", "data.frame"))
    rownames(obs) <- paste(modalityname, rownames(obs), sep="-")

  colnames(X) <- rownames(obs)
  rownames(X) <- rownames(var)

  raw <- NULL
  if ("raw" %in% names(root)) {
    raw <- root[['raw']]
    raw.X <- read_matrix(raw[['X']])
    raw.var <- read_table(raw[['var']])
    rownames(raw.X) <- rownames(raw.var)
    colnames(raw.X) <- colnames(X)
    if (nrow(raw.X) != nrow(X)) {
      warning(paste0("Only a subset of mod/", modalityname, "/raw/X is loaded, variables (features) that are not present in mod/", modalityname, "/X are discarded."))
      raw.X <- raw.X[rownames(X),]
    }
  }

  layers <- NULL
  custom_layers <- NULL
  if ("layers" %in% names(root)) {
    layers <- lapply(root[['layers']]$names, function(layer_name) {
      layer <- read_matrix(root[['layers']][[layer_name]])
      rownames(layer) <- rownames(X)
      colnames(layer) <- colnames(X)
      layer
    })
    names(layers) <- root[['layers']]$names
    custom_layers <- names(layers)[!names(layers) %in% c("counts")]
    if (length(custom_layers) > 0) {
      missing_on_read(paste0("some of mod/", modalityname, "/layers"), "custom layers, unless labeled 'counts'")
    }
  }

  # Assumptions:
  #   1. X -> counts
  #   2. raw & X -> data & scale.data
  #   3. layers['counts'] & X -> counts & data
  #   4. layers['counts'], raw, X -> counts, data, scale.data
  counts_as_layer <- !is.null(layers) && "counts" %in% names(layers)
  if (!counts_as_layer && is.null(raw)) {
    # 1
    assay <- Seurat::CreateAssayObject(counts = X)
  } else {
    if (!is.null(raw)) {
      if (counts_as_layer) {
        # 4
        assay <- Seurat::CreateAssayObject(counts = layers[['counts']])
        assay@data <- raw.X
        assay@scale.data <- X
      } else {
        # 2
        assay <- Seurat::CreateAssayObject(data = raw.X)
        assay@scale.data <- X
      }
    } else {
      # 3
      assay <- Seurat::CreateAssayObject(counts = layers[['counts']])
      assay@data <- X
    }
  }

  assay
}

read_attr_m <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_table(root[[attr_name]])
    dim_names <- rownames(attr_df)
  }
  attrm_name <- paste0(attr_name, "m")

  attrm <- list()
  if (attrm_name %in% names(root)) {
    attrm <- lapply(names(root[[attrm_name]]), function(space) {
      mx <- t(root[[attrm_name]][[space]]$read())
      if (dim(mx)[1] == 1) {
        mx <- t(mx)
      }
      rownames(mx) <- dim_names
      mx
    })

    names(attrm) <- names(root[[attrm_name]])
  }

  attrm
}

#' @import Seurat methods
read_attr_p <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_table(root[[attr_name]])
    dim_names <- rownames(attr_df)
  }
  attrp_name <- paste0(attr_name, "p")

  attrp <- list()
  if (attrp_name %in% names(root)) {
    attrp <- lapply(names(root[[attrp_name]]), function(graph) {
      mx <- read_matrix(root[[attrp_name]][[graph]])
      rownames(mx) <- dim_names
      colnames(mx) <- dim_names
      # Prevent automatic coersion based on equal dimensions
      if ("dsCMatrix" %in% class(mx)) {
        mx <- as(mx, "dgCMatrix")
      }
    })

    names(attrp) <- names(root[[attrp_name]])
  }

  attrp
}

# For some common reductions,
# there are conventional names for the loadings slots
OBSM2VARM <- list("X_pca" = "PCs", "X_mofa" = "LFs")
