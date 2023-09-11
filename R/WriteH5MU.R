#' @rdname WriteH5MU
setGeneric("WriteH5MU", function(object, file, overwrite = TRUE) standardGeneric("WriteH5MU"))

#' @rdname WriteH5AD
setGeneric("WriteH5AD", function(object, file, assay = NULL, overwrite = TRUE) standardGeneric("WriteH5AD"))

#' A helper function to write a modality (an assay) to an .h5mu file
#'
#' @keywords internal
#'
#' @import hdf5r methods
#' @importFrom Matrix t
WriteH5ADHelper <- function(object, assay, root, global = FALSE) {
  mod_object <- Seurat::GetAssay(object, assay)

  # .obs
  obs_names <- colnames(object)
  # There is no local metadata in Seurat objects
  if (global) {
    obs <- object@meta.data
  } else {
    obs <- data.frame(row.names = obs_names)
  }
  write_data_frame(root, "obs", obs)

  # .var
  var <- mod_object@meta.features
  var_names <- rownames(mod_object@meta.features)

  # Define highly variable features, if any
  if ('var.features' %in% slotNames(mod_object)) {
    if (length(mod_object@var.features) > 0) {
      var$highly_variable <- rownames(var) %in% mod_object@var.features
      message("Added .var['highly_variable'] with highly variable features")
    }
  }

  write_data_frame(root, "var", var)

  # .X, .layers['counts']. .raw.X
  # Assumptions:
  #   1. counts/data/scale.data -> X
  #   3. counts & data -> layers['counts'], X
  #   2. data & scale.data -> layers['data'], X
  #   4. counts & scale.data -> layers['counts'], X
  #   5. counts & data & scale.data -> layers['counts'], layers['data'], X

  x_names <- list("counts", "data", "scale.data")

  x <- lapply(x_names, function(x_name) {
    x <- NULL
    if (x_name %in% slotNames(mod_object)) {
      x <- Seurat::GetAssayData(mod_object, x_name)
      if (nrow(x) == 0 || ncol(x) == 0)
        x <- NULL
    }
    x
  })
  names(x) <- x_names

  if (!any(vapply(x, is.null, TRUE))) {
    # 5
    layers_group <- root$create_group("layers")
    write_matrix(layers_group, "counts", x[["counts"]])
    write_matrix(layers_group, "data", x[["data"]])
    write_matrix(root, "X", x[["scale.data"]])
  } else if (!is.null(x[["counts"]]) && !is.null(x[["scale.data"]])) {
    # 4
    layers_group <- root$create_group("layers")
    write_matrix(layers_group, "counts", x[["counts"]])
    write_matrix(root, "X", x[["scale.data"]])
  } else if (!is.null(x[["data"]]) && !is.null(x[["scale.data"]])) {
    # 3
    layers_group <- root$create_group("layers")
    write_matrix(layers_group, "data", x[["data"]])
    write_matrix(root, "X", x[["scale.data"]])
  } else if (!is.null(x[["counts"]]) && !is.null(x[["data"]])) {
    # 2
    layers_group <- root$create_group("layers")
    write_matrix(layers_group, "counts", x[["counts"]])
    write_matrix(root, "X", x[["data"]])
  } else {
    which_x <- which(!is.null(x))
    write_matrix(root, "X", x[[which_x]])
  }

  uns_group <- root$create_group("uns")

  # reductions -> .obsm
  if ('reductions' %in% slotNames(object)) {
    obsm_group <- root$create_group("obsm")
    varm_group <- root$create_group("varm")
    for (red_name in names(object@reductions)) {
      red <- object@reductions[[red_name]]
      emb <- t(red@cell.embeddings)
      emb_assay <- red@assay.used
      loadings <- red@feature.loadings

      modality_specific <- FALSE
      # Modality-specific reductions can be identified with all their feature names
      # coming from the @assay.used.
      if (!modality_specific) {
        if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
          if (all(rownames(loadings) %in% var_names)) {
            modality_specific <- TRUE
          }
        }
      }

      # Modality-specific reductions in Seurat objects
      # can also start with modality name by the convention used in this package.
      # Multimodal reductions also have the @assay.used set because this is enforced
      # by the current Seurat package.
      if (!is.null(emb_assay) && emb_assay != "" && emb_assay == assay) {
        # Only count reduction as modality-specific
        # if its name can be found in the reduction name or reduction key.
        # This is required since Seurat does require having an existing modality
        # in assay.used, which complicates loading multimodal embeddings.
        # The latter are currently loaded with the default assay set as assay.used.
        if (grepl(tolower(emb_assay), tolower(red_name)) || grepl(tolower(emb_assay), tolower(red@key))) {
          modality_specific <- TRUE
        }

        # Strip away modality name if the embedding starts with it
        if (emb_assay == substr(red_name, 1, length(emb_assay))) {
          red_name <- substr(red_name, length(emb_assay) + 1, length(red_name))
        }
      }

      if (!modality_specific) {
        next
      }

      write_matrix(obsm_group, paste0("X_", red_name), emb)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key = OBSM2VARM[[paste0("X_", red_name)]]
        }

        # If only a subset of features was used,
        # this has to be accounted for
        if (nrow(loadings) < nrow(var)) {
          warning(paste0("Loadings for ", red_name, " are computed only for a some features.",
            " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data."))
          all_loadings <- matrix(
            ncol = ncol(loadings),
            nrow = nrow(var)
          )
          rownames(all_loadings) <- rownames(var)
          all_loadings[rownames(loadings),] <- loadings
        } else {
          all_loadings <- loadings
        }

        write_matrix(varm_group, varm_key, t(all_loadings))
      }

      # stdev -> .uns[...]['variance']
      if (length(red@stdev) > 0) {
        if (!red_name %in% names(uns_group)) {
          uns_red <- uns_group$create_group(red_name)
          write_attribute(uns_red, "encoding-type", "dict")
          write_attribute(uns_red, "encoding-version", "0.1.0")
          write_matrix(uns_red, "variance", red@stdev^2)
        }
      }
    }
  }

  # graphs -> .obsp
  if ('graphs' %in% slotNames(object)) {
    obsp_group <- root$create_group("obsp")
    for (graph_name in names(object@graphs)) {
      graph <- object@graphs[[graph_name]]
      # Only write the graphs with the correct assay.used
      if ('assay.used' %in% slotNames(graph)) {
        if (length(graph@assay.used) > 0 && graph@assay.used == assay) {
          # Strip away modality name if the graph name starts with it:
          # RNA_distances -> distances
          if (assay == substr(graph_name, 1, nchar(graph@assay.used))) {
            graph_name <- substr(graph_name, nchar(graph@assay.used) + 1, nchar(graph_name))
            # Account for _, which is added by ReadH5AD / ReadH5MU
            if (substr(graph_name, 1, 1) == "_") {
              graph_name <- substr(graph_name, 2, nchar(graph_name))
            }
          }
          write_matrix(obsp_group, graph_name, graph)
        }
      }
    }
  }

  finalize_anndata_internal(root)

  TRUE
}

#' Write one assay to .h5ad
#'
#' This function writes the data of one of the assays (modalities) of a \code{Seurat} object into an .h5ad file.
#' The behavior of this function if NAs are present is undefined.
#'
#' The following slots are saved: count matrices (`@counts`, `@scale.data` and `@data`), `@metadata`, `@reductions`, `@feature.loadings`, `@graphs`.
#'
#' @param object \code{Seurat} object.
#' @param file Path to the .h5ad file.
#' @param assay Assay to write; can be omitted if there is a single assay in the object.
#' @param overwrite Boolean value to indicate if to overwrite the \code{file} if it exists (\code{TRUE} by default).
#'
#' @rdname WriteH5AD
#'
#' @import hdf5r
#'
#' @exportMethod WriteH5AD
setMethod("WriteH5AD", "Seurat", function(object, file, assay = NULL, overwrite = TRUE) {
  if (isFALSE(overwrite) && file.exists(file)) {
    stop(paste0("File ", file, " already exists. Use `overwrite = TRUE` to overwrite it or choose a different file name."))
  }

  h5 <- open_h5(file)

  # When multiple modalities are present,
  # an assay has to be specified.
  # Do not default to Seurat::DefaultAssay(object)
  # as it is not explicit, is hard to reason about,
  # and does not mean anything for MuData.
  if (length(object@assays) > 1 && is.null(assay)) {
    h5$close()
    stop(paste0(
      "An assay to be written has to be provided, one of: ",
      paste(names(object@assays), collapse = ", "),
      ".\nUse WriteH5MU() to write all the modalities."
    ))
  } else {
    assay <- names(object@assays)[1]
  }

  # "Global" attributes such as metadata have to be written
  WriteH5ADHelper(object, assay, h5, global = TRUE)

  finalize_anndata(h5)

  invisible(TRUE)
})

#' Create an .h5mu file with data from a \code{\link{Seurat}} object
#'
#' Save \code{\link{Seurat}} object to .h5mu file.
#' The behavior of this function if NAs are present is undefined.
#'
#' The following slots are saved: count matrices (`@counts`, `@scale.data` and `@data`), `@metadata`, `@reductions`, `@feature.loadings`, `@graphs`.
#'
#' @param object \code{Seurat} object.
#' @param file Path to the .h5mu file.
#' @param overwrite Boolean value to indicate if to overwrite the \code{file} if it exists (\code{TRUE} by default).
#'
#' @rdname WriteH5MU
#'
#' @import hdf5r methods
#'
#' @exportMethod WriteH5MU
setMethod("WriteH5MU", "Seurat", function(object, file, overwrite) {
  h5 <- open_h5(file)

  # .obs
  obs <- object@meta.data

  write_data_frame(h5, "obs", obs)

  modalities <- Seurat::Assays(object)

  h5mod <- h5$create_group("mod")
  h5mod$create_attr("mod-order", modalities)
  var_names <- lapply(modalities, function(mod) {
    mod_group <- h5$create_group(paste0("mod/", mod))

    WriteH5ADHelper(object, mod, mod_group)

    mod_object <- object[[mod]]
    rownames(mod_object)
  })
  names(var_names) <- modalities
  write_data_frame(h5, "var", do.call(c, var_names))

  uns_group <- h5$create_group("uns")

  # reductions -> .obsm
  # Reductions starting with modality name
  # that corresponds to the assay.used value
  # will be stored in .obsm slots of individual modalities:
  # RNAUMAP -> /mod/RNA/obsm/UMAP
  if ('reductions' %in% slotNames(object)) {
    for (red_name in names(object@reductions)) {
      red <- object@reductions[[red_name]]
      emb <- t(red@cell.embeddings)
      assay_emb <- red@assay.used
      loadings <- red@feature.loadings

      modality_specific <- FALSE
      # Modality-specific reductions can be identified with all their feature names
      # coming from the @assay.used.
      if (!modality_specific) {
        if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
          if (all(rownames(loadings) %in% var_names[[assay_emb]])) {
            modality_specific <- TRUE
          }
        }
      }

      # Modality-specific reductions in Seurat objects
      # can also start with modality name by the convention used in this package.
      # Multimodal reductions also have the @assay.used set because this is enforced
      # by the current Seurat package.
      if (!is.null(assay_emb) && assay_emb != "" && assay_emb %in% modalities) {
        # Only count reduction as modality-specific
        # if its name can be found in the reduction name or reduction key.
        # This is required since Seurat does require having an existing modality
        # in assay.used, which complicates loading multimodal embeddings.
        # The latter are currently loaded with the default assay set as assay.used.
        if (grepl(tolower(assay_emb), tolower(red_name)) || grepl(tolower(assay_emb), tolower(red@key))) {
          modality_specific <- TRUE
        }

        # Strip away modality name if the embedding starts with it
        if (assay_emb == substr(red_name, 1, length(assay_emb))) {
          red_name <- substr(red_name, length(assay_emb) + 1, length(red_name))
        }
      }

      if (modality_specific) {
        next
      }

      if (!"obsm" %in% names(h5)) {
        obsm <- h5$create_group("obsm")
      } else {
        obsm <- h5[["obsm"]]
      }

      write_matrix(obsm, paste0("X_", red_name), emb)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key <- OBSM2VARM[[paste0("X_", red_name)]]
        }

        if (modality_specific) {
          # this should have been written with WriteH5ADHelper
          next
        }

        if (!"varm" %in% names(h5)) {
          varm <- h5$create_group("varm")
        } else {
          varm <- h5[["varm"]]
        }

        # If only a subset of features was used,
        # this has to be accounted for
        var_names_for_loadings <- do.call(c, var_names)

        if (nrow(loadings) < length(var_names_for_loadings)) {
          warning(paste0("Loadings for ", red_name, " are computed only for a some features.",
            " For it, an array with full var dimension will be recorded as it has to be match the var dimension of the data."))
          all_loadings <- matrix(
            ncol = ncol(loadings),
            nrow = length(var_names_for_loadings)
          )
          rownames(all_loadings) <- var_names_for_loadings
          all_loadings[rownames(loadings),] <- loadings
        } else {
          all_loadings <- loadings
        }

        write_matrix(varm, varm_key, t(all_loadings))
      }

      # stdev -> .uns[...]['variance']
      if (length(red@stdev) > 0) {
        if (modality_specific) {
          # REMOVE: this should have been written with WriteH5ADHelper
          if (!red_name %in% names(h5[[paste0("mod/", assay_emb, "/uns")]])) {
            uns <- h5$create_group(paste0("mod/", assay_emb, "/uns/", red_name))
            write_attribute(uns, "encoding-type", "dict")
            write_attribute(uns, "encoding-version", "0.1.0")
          } else {
            uns <- uns_group[[paste0("mod/", assay_emb, "/uns/", red_name)]]
          }
        } else {
          if (!red_name %in% names(uns_group)) {
            uns <- uns_group$create_group(red_name)
            write_attribute(uns, "encoding-type", "dict")
            write_attribute(uns, "encoding-version", "0.1.0")
          } else {
            uns <- uns_group[[red_name]]
          }
        }
        write_matrix(uns, "variance", red@stdev^2)
      }
    }
  }

  # graphs -> .obsp
  if ('graphs' %in% slotNames(object)) {
    obsp_group <- h5$create_group("obsp")
    for (graph_name in names(object@graphs)) {
      graph <- object@graphs[[graph_name]]

      # Only write the graphs with no (correct) assay.used
      graph_no_assay <- FALSE
      if (!'assay.used' %in% slotNames(graph)) {
        graph_no_assay <- TRUE
      } else {
        if (!graph@assay.used %in% modalities)
          graph_no_assay <- TRUE
      }

      if (graph_no_assay) {
        write_matrix(obsp_group, graph_name, graph)
      }
    }
  }

  finalize_mudata(h5)

  invisible(TRUE)
})
