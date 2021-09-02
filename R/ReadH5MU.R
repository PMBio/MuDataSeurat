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
    # TODO: .raw
    # TODO: .layers

    var <- read_with_index(view[['var']])

    obs <- read_with_index(view[['obs']])
    if (is("obs", "data.frame"))
      rownames(obs) <- paste(mod, rownames(obs), sep="-")

    colnames(X) <- rownames(obs)
    rownames(X) <- rownames(var)

    assay <- Seurat::CreateAssayObject(counts = X)

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
