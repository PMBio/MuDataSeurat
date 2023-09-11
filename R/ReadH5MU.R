#' Read an .h5mu file and create a \code{\link{Seurat}} object.
#'
#' @param file Path to the .h5mu file.
#'
#' @return A \code{\link{Seurat}} object
#'
#' @export
ReadH5AD <- function(file) {
  # Connect to the the file
  h5 <- open_anndata(file)

  # Get metadata
  obs <- read_table(h5[["obs"]])
  var <- read_table(h5[["var"]])

  # X
  assay <- read_layers_to_assay(h5)

  # obsm
  obsm <- read_attr_m(h5, 'obs')

  # varm
  varm <- read_attr_m(h5, 'var')

  # obsp
  obsp <- read_attr_p(h5, 'obs')

  # If there are var pairs, there's no place to store it
  # in the Seurat object
  # var_pairs <- read_attr_p(h5, 'var')
  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # Create a Seurat object
  # If read from .h5mu modality, give an assay name
  path_fragments <- strsplit(file, "\\.h5mu")[[1]]
  if (length(path_fragments) == 2) {
    mod_path_fragments <- strsplit(path_fragments[2], "\\/")[[1]]
    assay_name <- mod_path_fragments[length(mod_path_fragments)]
    srt <- Seurat::CreateSeuratObject(assay, assay = assay_name)
  } else {
    srt <- Seurat::CreateSeuratObject(assay)
  }

  # Specify highly variable features
  if ("highly_variable" %in% colnames(var)) {
    if (is.logical(var$highly_variable)) {
      srt@assays[[1]]@var.features <- rownames(srt)[var$highly_variable]
    }
  }

  # Add metadata
  meta_data_names <- rownames(srt@meta.data)
  srt@meta.data <- cbind.data.frame(obs, srt@meta.data)
  rownames(srt@meta.data) <- meta_data_names

  # Add feature metadata
  meta_features_names <- rownames(srt@assays[[1]]@meta.features)
  srt@assays[[1]]@meta.features <- cbind.data.frame(var, srt@assays[[1]]@meta.features)
  rownames(srt@assays[[1]]@meta.features) <- meta_features_names

  # Add embeddings
  for (emb in names(obsm)) {
    emb_name <- gsub('X_', '', emb)

    maybe_loadings <- matrix()
    if (emb %in% names(OBSM2VARM)) {
      varm_key = OBSM2VARM[[emb]]
      if (varm_key %in% names(varm)) {
        maybe_loadings <- varm[[varm_key]]
      }
    }

    emb_stdev <- numeric()
    if ("uns" %in% names(h5)) {
      if (emb_name %in% names(h5[["uns"]])) {
        if ("variance" %in% names(h5[["uns"]][[emb_name]])) {
          emb_stdev <- sqrt(h5[["uns"]][[emb_name]][["variance"]]$read())
        }
      }
    }

    srt[[emb_name]] <- Seurat::CreateDimReducObject(
      embeddings = obsm[[emb]][rownames(obs),,drop=FALSE],
      loadings = maybe_loadings,
      key = paste0(emb_name, "_"),
      assay = Seurat::DefaultAssay(srt),
      stdev = emb_stdev
    )
  }

  # Add graphs
  srt@graphs <- lapply(obsp, Seurat::as.Graph)

  # Close the connection
  h5$close()

  srt
}

#' Create a \code{Seurat} object from .h5mu file contents
#'
#' @param file Path to the .h5mu file
#'
#' @import hdf5r Matrix Seurat
#' @importFrom utils hasName
#'
#' @return A \code{Seurat} object
#' '
#' @export ReadH5MU
ReadH5MU <- function(file) {
  # Connect to the the file
  h5 <- open_and_check_mudata(file)

  # Get assays (modalities)
  assays <- h5[["mod"]]$names
  if ("mod-order" %in% names(h5attributes(h5[["mod"]]))) {
    modorder <- h5attributes(h5[["mod"]])$`mod-order`
    if (all(assays %in% modorder)) {
      modorder <- modorder[modorder %in% assays]
      if (!any(duplicated(modorder))) {
        assays <- modorder
      }
    }
  }

  # Get global metadata
  metadata <- read_table(h5[["obs"]])

  # NOTE: there's no global feature metadata in the Seurat object
  ft_metadata <- tryCatch({
      read_table(h5[["var"]])
    },
    error = function(err) {
      warning(err)
      read_table(h5[["var"]], set_index = FALSE)
    }
  )

  if (ncol(ft_metadata) > 0)
    missing_on_read("/var", paste0("global variables metadata (", paste(colnames(ft_metadata), collapse = ", "), ")"))

  # Get (multimodal) embeddings
  embeddings <- read_attr_m(h5, 'obs', rownames(metadata))
  # If obs->mod mappings are in the file, dismiss them
  embeddings <- embeddings[!names(embeddings) %in% assays]

  # Get (multimodal) loadings
  # NOTE: features can be set as row names only if they are unique
  loadings <- read_attr_m(h5, 'var', rownames(ft_metadata))

  # Get obs pairs
  obs_pairs <- read_attr_p(h5, 'obs')

  # If there are var pairs, there's no place to store it
  # in the Seurat object
  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # mod/.../X, raw, and layers
  modalities <- lapply(assays, function(mod) {
    assay <- read_layers_to_assay(h5[['mod']][[mod]], mod)
  })
  names(modalities) <- assays

  # mod/.../obs
  mod_obs <- lapply(assays, function(mod) {
    read_table(h5[['mod']][[mod]][['obs']])
  })
  names(mod_obs) <- assays

  # mod/.../var
  mod_var <- lapply(assays, function(mod) {
    read_table(h5[['mod']][[mod]][['var']])
  })
  names(mod_var) <- assays

  # mod/.../obsm
  mod_obsm <- lapply(assays, function(mod) {
    read_attr_m(h5[['mod']][[mod]], 'obs')
  })
  names(mod_obsm) <- assays

  # mod/.../varm
  mod_varm <- lapply(assays, function(mod) {
    read_attr_m(h5[['mod']][[mod]], 'var')
  })
  names(mod_varm) <- assays

  # mod/.../obsp
  mod_obsp <- lapply(assays, function(mod) {
    read_attr_p(h5[['mod']][[mod]], 'obs')
  })
  names(mod_obsp) <- assays

  # If there are var pairs in individual modalities,
  # there's no place to store it in the Seurat object.
  for (mod in assays) {
    if ("varp" %in% names(h5[['mod']][[mod]])) {
      if (length(h5[['mod']][[mod]][['varp']]) > 0) {
        missing_on_read(paste0("/mod", mod, "/varp"), "pairwise annotation of variables")
      }
    }
  }

  var_pairs_names <- c()
  if ("varp" %in% names(h5))
    var_pairs_names <- names(h5[["varp"]])
  if (!is.null(var_pairs_names) && !isFALSE(var_pairs_names) && length(var_pairs_names) > 0)
    missing_on_read("/varp", "pairwise annotation of variables")

  # Only common observations can be read
  obs_names <- Reduce(intersect, lapply(modalities, colnames))
  mods_n_obs <- unique(vapply(modalities, ncol, 1))
  if (length(mods_n_obs) > 1 || mods_n_obs[1] != length(obs_names)) {
    warning("Only the intersection of observations (samples) is loaded. Observations that are not present in all the modalities (assays) are discarded.")
  }

  # Create a Seurat object
  srt <- Seurat::CreateSeuratObject(subset(modalities[[1]], cells = obs_names), assay = names(modalities)[1])
  for (modality in names(modalities)[2:length(modalities)]) {
    srt[[modality]] <- subset(modalities[[modality]], cells = obs_names)
  }

  # Metadata, features metadata, and variable features
  for (modality in names(modalities)) {
    # Append modality metadata
    srt@meta.data <- cbind.data.frame(srt@meta.data, mod_obs[[modality]][obs_names,])

    # Add modality feature metadata
    meta_features_names <- rownames(srt[[modality]]@meta.features)
    srt[[modality]]@meta.features <- cbind.data.frame(mod_var[[modality]], srt[[modality]]@meta.features)
    rownames(srt[[modality]]@meta.features) <- meta_features_names

    # Specify highly variable features
    if ("highly_variable" %in% colnames(mod_var[[modality]])) {
      if (is.logical(mod_var[[modality]]$highly_variable)) {
        srt[[modality]]@var.features <- rownames(srt[[modality]])[mod_var[[modality]]$highly_variable]
      }
    }
  }

  # Add joint embeddings
  for (emb in names(embeddings)) {
    emb_name <- toupper(gsub('X_', '', emb))

    maybe_loadings <- matrix()
    varm_key <- emb_name
    if (emb %in% names(OBSM2VARM)) {
      varm_key <- OBSM2VARM[[emb]]
    }
    if (hasName(loadings, varm_key)) {
      maybe_loadings <- loadings[[varm_key]]
    }

    emb_stdev <- numeric()
    if ("uns" %in% names(h5)) {
      if (emb_name %in% names(h5[["uns"]])) {
        if ("variance" %in% names(h5[["uns"]][[emb_name]])) {
          emb_stdev <- sqrt(h5[["uns"]][[emb_name]][["variance"]]$read())
        }
      }
    }

    srt[[emb_name]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings[[emb]][obs_names,,drop=FALSE],
      loadings = maybe_loadings,
      key = paste0(emb_name, "_"),
      stdev = emb_stdev,
      assay = Seurat::DefaultAssay(srt),  # this is not true but an existing assay must be provided
    )
  }

  # Do embeddings across modalities have unique names?
  # If each modality has e.g. X_pca, we will need to harmonise the names
  # by prepending modality name: e.g. RNAPCA.
  unique_emb <- FALSE
  all_embeddings <- c(names(embeddings), unlist(lapply(mod_obsm, function(obsm) names(obsm))))
  if (!any(duplicated(all_embeddings))) {
    unique_emb <- TRUE
  }

  # Add modality-specific embeddings
  for (mod in names(mod_obsm)) {
    mod_embeddings <- mod_obsm[[mod]]
    for (emb in names(mod_embeddings)) {
      emb_name <- gsub('X_', '', emb)
      if (!startsWith(emb_name, mod))
        emb_name <- toupper(emb_name)
      modemb_name <- emb_name
      if (!unique_emb)
        modemb_name <- paste(mod, emb_name, sep = "")
      # Embeddings keys will have to follow the format alphanumericcharacters_, e.g. RNAPCA_.

      maybe_loadings <- matrix()
      varm_key <- emb_name
      if (emb %in% names(OBSM2VARM))
        varm_key = OBSM2VARM[[emb]]
      if (hasName(mod_varm[[mod]], varm_key))
        maybe_loadings <- mod_varm[[mod]][[varm_key]]

      emb_stdev <- numeric()
      h5_mod <- h5[["mod"]][[mod]]
      if ("uns" %in% names(h5_mod)) {
        if (emb_name %in% names(h5_mod[["uns"]])) {
          if ("variance" %in% names(h5_mod[["uns"]][[emb_name]])) {
            emb_stdev <- sqrt(h5_mod[["uns"]][[emb_name]][["variance"]]$read())
          }
        }
      }

      srt[[modemb_name]] <- Seurat::CreateDimReducObject(
        embeddings = mod_embeddings[[emb]][obs_names,,drop=FALSE],
        loadings = maybe_loadings,
        key = paste0(modemb_name, "_"),
        stdev = emb_stdev,
        assay = mod,
      )
    }
  }

  # Add graphs

  # Only take into account common observations
  if (length(obs_pairs) > 0) {
    srt@graphs <- lapply(obs_pairs, function(graph) {
      graph[obs_names,obs_names,drop=FALSE]
    })
    names(srt@graphs) <- names(obs_pairs)
  }

  for (mod in names(mod_obsp)) {
    for (graph in names(mod_obsp[[mod]])) {
      graph_name <- graph
      # /mod/RNA/obsp/distances -> @graphs$RNA_distances
      if (graph %in% names(srt@graphs)) {
        graph_name <- paste(mod, graph, sep = "_")
      }
      srt@graphs[[graph_name]] <- Seurat::as.Graph(mod_obsp[[mod]][[graph]][obs_names, obs_names])
      srt@graphs[[graph_name]]@assay.used <- mod
    }
  }


  # Close the connection
  h5$close_all()

  srt
}
