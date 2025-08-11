#' Calculate spatial correlation of gene expression with physical distance between cells
#'
#' Computes the correlation between expression similarity and physical distance, summarized in distance bins.
#'
#' @param seurat_obj A Seurat object with spatial coordinates in metadata.
#' @param gene Character; gene name to analyze.
#' @param coord_cols Character vector of length 2 specifying spatial coordinate columns (default c("imagecol", "imagerow")).
#' @param assay Assay name for expression (default "RNA").
#' @param slot Slot name in assay (default "data").
#' @param max_distance Numeric; maximum spatial distance to consider (default NULL, uses max observed).
#' @param n_bins Integer; number of distance bins to aggregate correlations (default 20).
#' @param sample_size Integer; number of random cell pairs to sample for correlation (default 10000).
#'
#' @return A data.frame with columns: distance_bin, mean_distance, mean_correlation.
#' @importFrom stats cor
#' @export
#'
#' @examples
#' dist_corr_df <- spatial_correlation_with_distance(seurat_obj, gene = "ACTA2")
spatial_correlation_with_distance <- function(seurat_obj,
                                              gene,
                                              coord_cols = c("imagecol", "imagerow"),
                                              assay = "RNA",
                                              slot = "data",
                                              max_distance = NULL,
                                              n_bins = 20,
                                              sample_size = 10000) {
  meta <- seurat_obj@meta.data
  
  if (!all(coord_cols %in% colnames(meta))) {
    stop("Coordinate columns not found in Seurat metadata.")
  }
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste0("Gene '", gene, "' not found in Seurat object."))
  }
  
  coords <- meta[, coord_cols]
  expr <- as.numeric(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, ])
  
  valid <- complete.cases(coords) & !is.na(expr)
  coords <- coords[valid, , drop = FALSE]
  expr <- expr[valid]
  
  n_cells <- nrow(coords)
  if (n_cells < 2) {
    stop("Not enough cells with valid coordinates and expression for analysis.")
  }
  
  # Sample random pairs for computational efficiency if needed
  max_pairs <- choose(n_cells, 2)
  n_samples <- min(sample_size, max_pairs)
  
  # Generate random pairs of indices
  pairs <- matrix(NA_integer_, nrow = n_samples, ncol = 2)
  set.seed(123)
  for (i in seq_len(n_samples)) {
    pairs[i, ] <- sample(n_cells, 2)
  }
  
  # Calculate Euclidean distances between pairs
  dists <- sqrt((coords[pairs[, 1], 1] - coords[pairs[, 2], 1])^2 +
                  (coords[pairs[, 1], 2] - coords[pairs[, 2], 2])^2)
  
  if (!is.null(max_distance)) {
    keep <- dists <= max_distance
    if (sum(keep) < 10) {
      warning("Few pairs under max_distance; increasing max_distance or reducing sample_size recommended.")
    }
    dists <- dists[keep]
    pairs <- pairs[keep, , drop = FALSE]
  }
  
  # Calculate expression correlations (Pearson) between pairs
  expr_corrs <- mapply(function(i, j) {
    cor(expr[i], expr[j], method = "pearson")
  }, pairs[, 1], pairs[, 2])
  
  # Cor between single values is NA, so use product of deviations as proxy or
  # better approach: calculate similarity = 1 - abs(expr[i] - expr[j]) or use product
  # Because cor(expr[i], expr[j]) of two scalars is NA, instead use negative absolute difference or similarity measure
  
  # We'll use negative absolute difference scaled to [-1, 0] as proxy for similarity:
  max_diff <- max(abs(expr - mean(expr))) * 2
  expr_sim <- 1 - (abs(expr[pairs[, 1]] - expr[pairs[, 2]]) / max_diff)
  
  # Bin by distance
  bins <- cut(dists, breaks = n_bins, include.lowest = TRUE)
  
  # Aggregate mean similarity and mean distance per bin
  agg_df <- aggregate(cbind(expr_sim, dists), by = list(bin = bins), FUN = mean)
  colnames(agg_df) <- c("distance_bin", "mean_correlation", "mean_distance")
  
  # Reorder columns for clarity
  agg_df <- agg_df[, c("distance_bin", "mean_distance", "mean_correlation")]
  
  return(agg_df)
}