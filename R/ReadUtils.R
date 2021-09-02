#' @import hdf5r
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

missing_on_read <- function(loc, desc = "") {
  details <- ""
  if (!is.null(desc) && desc != "") {
    details <- paste0("Seurat does not support ", desc, ".")
  }
  warning(paste0("Missing on read: ", loc, details))
}


read_with_index <- function(dataset) {
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
          values <- factor(as.integer(values), labels = values_labels$read()[1:length(unique(values))])
        }
      }
      values
    })
    table <- data.frame(Reduce(cbind, col_list))
    colnames(table) <- columns

    if (indexcol %in% colnames(table)) {
      rownames(table) <- table[,indexcol,drop=TRUE]
      table <- table[,!colnames(table) %in% c(indexcol),drop=FALSE]
    }

    # DEPRECATED:
    # For consistency with other tools, this is done via references (see above)
    # Make factors out of categorical data
    # if ("__categories" %in% names(dataset)) {
    #   cats <- dataset[["__categories"]]
    #   for (cat in names(cats)) {
    #     table[[cat]] <- factor(as.integer(table[[cat]]) + 1, labels = cats[[cat]]$read())
    #   }
    # }

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

    if (indexcol %in% colnames(table)) {
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
      if ("shape" %in% h5attr_names(dataset)) {
        X_dims <- h5attr(dataset, "shape")
      } else {
        X_dims <- c(max(i), max(p))
      }
      X <- Matrix::Matrix(0, X_dims[1], X_dims[2])
      X@i <- i
      X@p <- p
      X@x <- x
      t(X)
    } else {
      dataset$read()
    }
}

read_attr_m <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_with_index(root[[attr_name]])
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

read_attr_p <- function(root, attr_name, dim_names = NULL) {
  if (is.null(dim_names)) {
    attr_df <- read_with_index(root[[attr_name]])
    dim_names <- rownames(attr_df)
  }
  attrp_name <- paste0(attr_name, "p")

  attrp <- list()
  if (attrp_name %in% names(root)) {
    attrp <- lapply(names(root[[attrp_name]]), function(graph) {
      mx <- read_matrix(root[[attrp_name]][[graph]])
      rownames(mx) <- dim_names
      colnames(mx) <- dim_names
      mx
    })
    
    names(attrp) <- names(root[[attrp_name]])
  }

  attrp
}


OBSM2VARM <- list("X_pca" = "PCs", "X_mofa" = "LFs")
