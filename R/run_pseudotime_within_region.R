#' Run pseudotime inference separately within each spatial region
#'
#' Performs trajectory inference (pseudotime) within each spatial contour/region using
#' reduced dimensional embeddings or expression data.
#'
#' @param seurat_obj A Seurat object.
#' @param region_col Metadata column name defining spatial regions/contours (default "assigned_region").
#' @param reduction Character; name of dimensionality reduction to use (e.g., "pca", "umap") (default "pca").
#' @param dims Integer vector; dimensions of the reduction to use (default 1:10).
#' @param cluster_col Metadata column with cell cluster labels (optional; if NULL, no clustering).
#' @param pseudotime_col_prefix Prefix for new pseudotime metadata columns per region (default "pseudotime").
#' @param slingshot_args List of additional arguments to pass to `slingshot::slingshot()`.
#'
#' @return Updated Seurat object with added pseudotime metadata columns named `<pseudotime_col_prefix>_<region>`.
#' @importFrom slingshot slingshot
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
#'
#' @examples
#' seurat_obj <- run_pseudotime_within_region(seurat_obj, region_col = "assigned_region", reduction = "pca", dims = 1:5)
run_pseudotime_within_region <- function(seurat_obj,
                                         region_col = "assigned_region",
                                         reduction = "pca",
                                         dims = 1:10,
                                         cluster_col = NULL,
                                         pseudotime_col_prefix = "pseudotime",
                                         slingshot_args = list()) {
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("Package 'slingshot' is required but not installed.")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed.")
  }
  
  meta <- seurat_obj@meta.data
  if (!(region_col %in% colnames(meta))) {
    stop(paste0("Metadata column '", region_col, "' not found in Seurat object."))
  }
  
  if (!(reduction %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", reduction, "' not found in Seurat object."))
  }
  
  # Prepare new metadata list to store pseudotime columns
  new_meta <- meta
  
  # For each region, run pseudotime inference separately
  regions <- unique(meta[[region_col]])
  for (region in regions) {
    cells_in_region <- rownames(meta)[meta[[region_col]] == region]
    
    if (length(cells_in_region) < 10) {
      warning(paste("Region", region, "has fewer than 10 cells; skipping pseudotime inference."))
      next
    }
    
    # Extract embeddings for these cells and selected dims
    emb <- Embeddings(seurat_obj, reduction = reduction)[cells_in_region, dims, drop = FALSE]
    
    # Optional cluster labels for slingshot
    cluster_labels <- NULL
    if (!is.null(cluster_col)) {
      if (!(cluster_col %in% colnames(meta))) {
        warning(paste0("Cluster column '", cluster_col, "' not found; skipping clustering for region ", region))
      } else {
        cluster_labels <- meta[cells_in_region, cluster_col]
        if (length(unique(cluster_labels)) < 2) {
          warning(paste("Not enough clusters in region", region, "for slingshot; skipping clustering."))
          cluster_labels <- NULL
        }
      }
    }
    
    # Construct SingleCellExperiment object
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(emb)))
    # Set reducedDims slot (slingshot looks for 'reducedDims')
    SingleCellExperiment::reducedDims(sce) <- list(RED = emb)
    
    # Run slingshot
    args <- c(list(sce = sce, reducedDim = "RED"), slingshot_args)
    if (!is.null(cluster_labels)) {
      args$cluster <- cluster_labels
    }
    sce <- do.call(slingshot::slingshot, args)
    
    # Extract pseudotime: if multiple lineages, take the first one
    pseudotime <- slingshot::slingPseudotime(sce)
    if (is.matrix(pseudotime)) {
      pseudotime_vec <- pseudotime[, 1]
    } else {
      pseudotime_vec <- pseudotime
    }
    
    # Add pseudotime to metadata with unique column name per region
    col_name <- paste0(pseudotime_col_prefix, "_", region)
    new_meta[cells_in_region, col_name] <- pseudotime_vec
  }
  
  seurat_obj@meta.data <- new_meta
  
  return(seurat_obj)
}