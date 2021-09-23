.mudataversion <- "0.1.0"
.anndataversion <- "0.1.0"
.name <- paste0(getPackageName(), ".r")
.version <- as.character(packageVersion(getPackageName()))

#' @import hdf5r
open_h5 <- function(filename) {
    h5p_create <- H5P_FILE_CREATE$new()
    h5p_create$set_userblock(512)
    H5File$new(filename, mode="w", file_create_pl=h5p_create)
}

#' @import hdf5r
finalize_mudata <- function(h5) {
    h5$create_attr("encoding-type", "MuData", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .mudataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))

    filename <- h5$get_filename()
    h5$close_all()
    h5 <- file(filename, "r+b")
    writeChar(paste0("MuData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
    close(h5)
}

#' @import hdf5r
finalize_anndata_internal <- function(h5) {
    h5$create_attr("encoding-type", "AnnData", space=H5S$new("scalar"))
    h5$create_attr("encoding-version", .anndataversion, space=H5S$new("scalar"))
    h5$create_attr("encoder", .name, space=H5S$new("scalar"))
    h5$create_attr("encoder-version", .version, space=H5S$new("scalar"))
}

#' @import hdf5r
finalize_anndata <- function(h5, internal = FALSE) {
    if (internal) {
      finalize_anndata_internal(h5)
    }
    filename <- h5$get_filename()
    h5$close_all()
    h5 <- file(filename, "r+b")
    writeChar(paste0("AnnData (format-version=", .mudataversion, ";creator=", .name, ";creator-version=", .version, ")"), h5)
    close(h5)
}


write_data_frame <- function(attr_group, attr_df) {
  attr_columns <- colnames(attr_df)

  attr_df["_index"] <- rownames(attr_df)

  attr_df <- attr_df[,c("_index", attr_columns),drop=FALSE]

  categories <- list()
  for (col in colnames(attr_df)) {
    v <- attr_df[[col]]
    if ("factor" %in% class(v)) {
      # Write a factor
      categories[[col]] <- levels(v)
      attr_group$create_dataset(col, as.integer(v) - 1, dtype = h5types$H5T_NATIVE_INT)
    } else {
      attr_group$create_dataset(col, v)
    }
  }
  if (length(categories) > 0) {
    cats <- attr_group$create_group("__categories")
    for (cat in names(categories)) {
      cat_dataset <- cats$create_dataset(cat, categories[[cat]])
      cat_dataset$create_attr("ordered", FALSE, space = H5S$new("scalar"))
      attr_group[[cat]]$create_attr("categories", 
                                    cats$create_reference(cat), 
                                    space = H5S$new("scalar"))
    }
  }

  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-version", "0.1.0", space = H5S$new("scalar"))
  if (length(attr_columns) > 0) {
    attr_group$create_attr("column-order", attr_columns)
  } else {
    # When there are no columns, null buffer can't be written to a file.
    attr_group$create_attr("column-order", dtype=h5types$H5T_NATIVE_DOUBLE, space=H5S$new("simple", 0, 0))
  }
  
}

# Only write _index (obs_names or var_names)
write_names <- function(attr_group, attr_names) {
  attr_group$create_dataset("_index", attr_names)

  # Write attributes
  attr_group$create_attr("_index", "_index", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-type", "dataframe", space = H5S$new("scalar"))
  attr_group$create_attr("encoding-version", "0.1.0", space = H5S$new("scalar"))
  # When there are no columns, null buffer can't be written to a file.
  attr_group$create_attr("column-order", dtype=h5types$H5T_NATIVE_DOUBLE, space=H5S$new("simple", 0, 0))
  
}

write_sparse_matrix <- function(root, x, sparse_type) {
  root$create_dataset("indices", x@i)
  root$create_dataset("indptr", x@p)
  root$create_dataset("data", x@x)
  h5attr(root, "shape") <- dim(x)
  root$create_attr("encoding-type", sparse_type, space=H5S$new("scalar"))
  root$create_attr("encoding-version", "0.1.0", space=H5S$new("scalar"))
}
