#' @description Create a \code{Seurat} object from .h5mu file contents
#'
#' @import hdf5r, Matrix, Seurat
#'
#' @return A \code{Seurat} object with the data requested
#' '
#' @exportMethod ReadH5MU
ReadH5MU <- function(file) {
  # For some common reductions, 
  # there are conventional names for the loadings slots
  OBSM2VARM <- list("X_pca" = "PCs", "X_mofa" = "LFs")
    
  # Connect to the the file
  h5 <- open_and_check_mudata(file)

  # Get assays (modalities)
  assays <- h5[['mod']]$names

  # Get global metadata
  metadata <- read_with_index(h5[["obs"]])
  ft_metadata <- read_with_index(h5[["var"]])

  # Get (multimoal) embeddings
  embeddings <- read_attr_m(h5, 'obs', rownames(metadata))
  # If obs->mod mappings are in the file, dismiss them
  embeddings <- embeddings[!names(embeddings) %in% assays]

  # Get (multimoal) loadings
  loadings <- read_attr_m(h5, 'var', rownames(ft_metadata))

  # Get obs pairs
  # FIXME
  # obs_pairs <- read_attr_p(h5, 'obs')

  # If there are var pairs, there's no place to store it 
  # in the Seurat object
  var_pairs <- read_attr_p(h5, 'var')
  if (!is.null(var_pairs) && length(var_pairs) > 0) 
    missing_on_read("/varp", "pairwise annotation of variables")
  
  # mod/.../X
  modalities <- lapply(assays, function(mod) {
    view <- h5[['mod']][[mod]]

    X <- read_matrix(view[['X']])

    var <- read_with_index(view[['var']])

    obs <- read_with_index(view[['obs']])
    if (is("obs", "data.frame"))
      rownames(obs) <- paste(mod, rownames(obs), sep="-")

    colnames(X) <- rownames(obs)
    rownames(X) <- rownames(var)

    raw <- NULL
    if ("raw" %in% names(view)) {
      raw <- view[['raw']]
      raw.X <- read_matrix(raw[['X']])
      raw.var <- read_with_index(raw[['var']])
      rownames(raw.X) <- rownames(raw.var)
      colnames(raw.X) <- colnames(X)
      if (nrow(raw.X) != nrow(X)) {
        warning(paste0("Only a subset of mod/", mod, "/raw/X is loaded, variables (features) that are not present in mod/", mod, "/X are discarded."))
        raw.X <- raw.X[rownames(X),]
      }
    }

    layers <- NULL
    custom_layers <- NULL
    if ("layers" %in% names(view)) {
      layers <- lapply(view[['layers']]$names, function(layer_name) {
        layer <- read_matrix(view[['layers']][[layer_name]])
        rownames(layer) <- rownames(X)
        colnames(layer) <- colnames(X)
        layer
      })
      names(layers) <- view[['layers']]$names
      custom_layers <- names(layers)[!names(layers) %in% c("counts")]
      if (length(custom_layers) > 0) {
        missing_on_read(paste0("some of mod/", mod, "/layers"), "custom layers, unless labeled 'counts'")
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
          # 2
          assay <- Seurat::CreateAssayObject(data = raw.X)
          assay@scale.data <- X 
        } else {
          # 4
          assay <- Seurat::CreateAssayObject(counts = layers[['counts']])
          assay@data <- raw.X
          assay@scale.data <- X 
        }
      } else {
        # 3
        assay <- Seurat::CreateAssayObject(counts = layers[['counts']])
        assay@data <- X
      }
    }

    assay
  })
  names(modalities) <- assays

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

  # TODO: .obsp
  # TODO: .varp

  # Only common observations can be read
  obs_names <- Reduce(intersect, lapply(modalities, colnames))
  mods_n_obs <- unique(vapply(modalities, ncol, 1))
  if (length(mods_n_obs) > 1 || length(mods_n_obs[1]) != length(obs_names)) {
    warning("Only the intersection of observations (samples) is loaded. Observations that are not present in all the modalities (assays) are discarded.")
  }

  # Create a Seurat object
  srt <- Seurat::CreateSeuratObject(modalities[[1]][,obs_names], assay = names(modalities)[1])
  for (modality in names(modalities)[2:length(modalities)]) {
    srt[[modality]] <- subset(modalities[[modality]], cells = obs_names)
  }

  # Add joint embeddings
  for (emb in names(embeddings)) {
    emb_name <- toupper(gsub('X_', '', emb))

    maybe_loadings <- matrix()
    if (emb %in% names(OBSM2VARM)) {
      varm_key = OBSM2VARM[[emb]]
      maybe_loadings <- loadings[[varm_key]]
    } 
    srt[[emb_name]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings[[emb]][obs_names,,drop=FALSE], 
      loadings = maybe_loadings,
      key = paste0(emb_name, "_"),
      assay = Seurat::DefaultAssay(srt),  # this is not true but an existing assay must be provided
    )
  }

  # Add modality-specific embeddings
  for (mod in names(mod_obsm)) {
    mod_embeddings <- mod_obsm[[mod]]
    for (emb in names(mod_embeddings)) {
      emb_name <- paste(mod, toupper(gsub('X_', '', emb)), sep = "")
      # Embeddings keys will have to follow the format alphanumericcharacters_, e.g. RNAPCA_.

      maybe_loadings <- matrix()
      if (emb %in% names(OBSM2VARM)) {
        varm_key = OBSM2VARM[[emb]]
        maybe_loadings <- mod_varm[[mod]][[varm_key]]
      }

      srt[[emb_name]] <- Seurat::CreateDimReducObject(
        embeddings = mod_embeddings[[emb]][obs_names,,drop=FALSE],
        loadings = maybe_loadings,
        key = paste0(emb_name, "_"), 
        assay = mod,
      )
    }
  }

  # Add graphs
  # FIXME
  # if (length(obs_pairs) > 0) {
  #   srt@graphs <- lapply(obs_pairs, function(graph) {
  #     graph[obs_names,obs_names,drop=FALSE]
  #   })
  #   names(srt@graphs) <- names(obs_pairs)
  # }

  # TODO: Data from .uns["neighbors"].

  # Close the connection
  h5$close_all()

  srt
}
