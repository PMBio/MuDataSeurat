setGeneric("WriteH5MU", function(object, file, overwrite = TRUE) standardGeneric("WriteH5MU"))
setGeneric("WriteH5AD", function(object, file, assay = NULL, overwrite = TRUE) standardGeneric("WriteH5AD"))

#' A helper function to write a modality (an assay) to an .h5mu file
#'
#' @keywords internal
#' 
#' @import hdf5r
#' @importFrom Matrix t
WriteH5ADHelper <- function(object, assay, root) {

  mod_object <- Seurat::GetAssay(object, assay)

  # .obs
  obs_group <- root$create_group("obs")
  # There is no local metadata in Seurat objects
  obs_names <- colnames(object)
  obs <- data.frame(row.names = obs_names)
  write_data_frame(obs_group, obs)

  # .var
  var <- mod_object@meta.features

  # Define highly variable features, if any
  if ('var.features' %in% slotNames(mod_object)) {
    if (length(mod_object@var.features) > 0) {
      message("Defining highly variable features...")
      var$highly_variable <- rownames(var) %in% mod_object@var.features
    }
  }

  var_group <- root$create_group("var")
  write_data_frame(var_group, var)

  # .X, .layers['counts']. .raw.X
  if ('counts' %in% slotNames(mod_object)) {
    x_counts <- Seurat::GetAssayData(mod_object, 'counts')
    if (nrow(x_counts) != 0 && ncol(x_counts) != 0) {
      sparse_type <- ifelse(class(x_counts) == "dgCMatrix", "csc_matrix", "csr_matrix")
      # case 1: only counts available
      if (!(('data' %in% slotNames(mod_object)) || ('scale.data' %in% slotNames(mod_object)))) {
        if ("i" %in% slotNames(x_counts)) {
          # sparse matrix
          x_counts <- Matrix::t(x_counts)
          counts_group <- root$create_group("X")
          write_sparse_matrix(counts_group, x_counts, sparse_type)
        } else {
          # dense matrix
          root$create_dataset("X", t(x_counts))
        }
      } else {
        layers_group <- root$create_group("layers")
        if ("i" %in% slotNames(x_counts)) {
          # sparse matrix
          x_counts <- Matrix::t(x_counts)
          counts_group <- layers_group$create_group("counts")
          write_sparse_matrix(counts_group, x_counts, sparse_type)
        } else {
          # dense matrix
          layers_group$create_dataset("counts", t(x_counts))
        }
        if ('data' %in% slotNames(mod_object)) {
          x_data <- Seurat::GetAssayData(mod_object, 'data')
          sparse_type <- ifelse(class(x_data) == "dgCMatrix", "csc_matrix", "csr_matrix")
          if ('scale.data' %in% slotNames(mod_object) && length(mod_object@scale.data) > 0) {
            # case 2: counts, data, and scale.data are available
            # .X
            x_scaled <- t(Seurat::GetAssayData(mod_object, 'scale.data'))
            root$create_dataset("X", x_scaled)
            # .raw
            raw_group <- root$create_group("raw")
            if ("i" %in% slotNames(x_data)) {
              # sparse matrix
              x_data <- Matrix::t(x_data)
              data_group <- raw_group$create_group("X")
              write_sparse_matrix(data_group, x_data, sparse_type)
            } else {
              # dense matrix
              raw_group$create_dataset("X", t(x_data))
            }
          } else {
            # case 3: counts and data are available but not scale.data
            if ("i" %in% slotNames(x_data)) {
              # sparse matrix
              x_data <- Matrix::t(x_data)
              data_group <- root$create_group("X")
              write_sparse_matrix(data_group, x_data, sparse_type)
            } else {
              # dense matrix
              root$create_dataset("X", t(x_data))
            }
          }
        }
        # 'data' should to be available when 'scale.data' is available
      }
    }
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

      if (!is.null(emb_assay) && emb_assay != "") {
        if (emb_assay != assay) {
          next
        }
        if (!is.null(loadings) && nrow(loadings) != nrow(object)) {
          next
        }

        # Strip away modality name if the embedding starts with it
        if (emb_assay == substr(red_name, 1, nchar(emb_assay))) {
          red_name <- substr(red_name, nchar(emb_assay) + 1, nchar(red_name))
        }
      }

      obsm_group$create_dataset(paste0("X_", red_name), emb)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key = OBSM2VARM[[paste0("X_", red_name)]]
        }

        varm_group$create_dataset(varm_key, t(loadings))
      }

      # stdev -> .uns[...]['variance']
      if (length(red@stdev) > 0) {
        if (!red_name %in% names(uns_group)) {
          uns_red <- uns_group$create_group(red_name)
          uns_red$create_dataset("variance", red@stdev ** 2)
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
          graph_group <- obsp_group$create_group(graph_name)
          write_sparse_matrix(graph_group, graph, "csc_matrix")
        }
      }
    }
  }

  finalize_anndata_internal(root)

  TRUE
}

#' Write one assay to .h5ad
#' 
#' This function writes the data of one of the assays (modalities) of a Seurat object 
#' into an AnnData-compatible .h5ad file.
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
  
  WriteH5ADHelper(object, assay, h5)
  
  finalize_anndata(h5)

  TRUE
})

#' Create an .h5mu file with data from a \code{\link{Seurat}} object
#' 
#' @description Save \code{\link{Seurat}} object to .h5mu file
#'
#' @import hdf5r
#'
#' @exportMethod WriteH5MU
setMethod("WriteH5MU", "Seurat", function(object, file, overwrite) {
  h5 <- open_h5(file)

  # .obs
  obs_group <- h5$create_group("obs")
  obs <- object@meta.data

  write_data_frame(obs_group, obs)

  modalities <- Seurat::Assays(object)

  h5$create_group("mod")
  var_names <- lapply(modalities, function(mod) {
    mod_group <- h5$create_group(paste0("mod/", mod))

    WriteH5ADHelper(object, mod, mod_group)

    mod_object <- object[[mod]]
    rownames(mod_object)
  })

  # global .var will only contain rownames
  # NOTE: creating a data.frame fails for objects 
  # that have the same feature name(s) across different modalities
  # var <- data.frame(row.names = do.call(c, var_names))
  # write_data_frame(var_group, var)
  var_group <- h5$create_group("var")
  write_names(var_group, do.call(c, var_names))

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
        # REMOVE: this should have been written with WriteH5ADHelper
        # Modality-specific obsm
        if (!"obsm" %in% names(h5[[paste0("mod/", assay_emb)]])) {
          obsm <- h5$create_group(paste0("mod/", assay_emb, "/obsm"))
        } else {
          obsm <- h5[[paste0("mod/", assay_emb, "/obsm")]]
        }
      } else {
        # Global obsm
        if (!"obsm" %in% names(h5)) {
          obsm <- h5$create_group("obsm")
        } else {
          obsm <- h5[["obsm"]]
        }
      }
      
      obsm$create_dataset(paste0("X_", red_name), emb)

      # loadings -> .varm
      if (!is.null(loadings) && ncol(loadings) == ncol(red)) {
        varm_key <- red_name
        if (paste0("X_", red_name) %in% names(OBSM2VARM)) {
          varm_key = OBSM2VARM[[red_name]]
        }

        if (modality_specific) {
          # REMOVE: this should have been written with WriteH5ADHelper
          if (!"varm" %in% names(h5[[paste0("mod/", assay_emb)]])) {
            varm <- h5$create_group(paste0("mod/", assay_emb, "/varm"))
          } else {
            varm <- h5[[paste0("mod/", assay_emb, "/varm")]]
          }
        } else {
          if (!"varm" %in% names(h5)) {
            varm <- h5$create_group("varm")
          } else {
            varm <- h5[["varm"]]
          }
        }
        varm$create_dataset(varm_key, t(loadings))
      }

      # stdev -> .uns[...]['variance']
      if (length(red@stdev) > 0) {
        if (modality_specific) {
          # REMOVE: this should have been written with WriteH5ADHelper
          if (!red_name %in% names(h5[[paste0("mod/", assay_emb, "/uns")]])) {
            uns <- h5$create_group(paste0("mod/", assay_emb, "/uns/", red_name))
          } else {
            uns <- uns_group[[paste0("mod/", assay_emb, "/uns/", red_name)]]
          }
        } else {
          if (!red_name %in% names(uns_group)) {
            uns <- uns_group$create_group(red_name)
          } else {
            uns <- uns_group[[red_name]]
          }
        }
        uns$create_dataset("variance", red@stdev ** 2)
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
        graph_group <- obsp_group$create_group(graph_name)
        write_sparse_matrix(graph_group, graph, "csc_matrix")
      }
    }
  }

  finalize_mudata(h5)

  TRUE
})
