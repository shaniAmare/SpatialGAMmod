#' Plot heatmaps of cell type proportions across regions
#'
#' Generates heatmaps showing proportions of cell types (or categories) across spatial regions/groups.
#'
#' @param prop_df Data frame with columns: region/group, cell type, and proportion.
#' @param region_col Name of column for regions/groups (default "assigned_region").
#' @param celltype_col Name of column for cell types (default "CellType").
#' @param proportion_col Name of numeric column with proportions (default "proportion").
#' @param cluster_rows Logical; whether to cluster heatmap rows (default TRUE).
#' @param cluster_cols Logical; whether to cluster heatmap columns (default TRUE).
#' @param scale Logical; whether to scale proportions by row or column ("row", "column", or FALSE) (default FALSE).
#' @param color_palette Color palette function or vector (default viridis).
#'
#' @return Invisibly returns the heatmap object (from pheatmap).
#' @importFrom pheatmap pheatmap
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' proportion_heatmaps(prop_df, region_col = "assigned_region", celltype_col = "CellType", proportion_col = "proportion")
proportion_heatmaps <- function(prop_df,
                                region_col = "assigned_region",
                                celltype_col = "CellType",
                                proportion_col = "proportion",
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                scale = FALSE,
                                color_palette = viridis::viridis) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required but not installed.")
  }
  
  # Check columns exist
  req_cols <- c(region_col, celltype_col, proportion_col)
  if (!all(req_cols %in% colnames(prop_df))) {
    stop("Data frame must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # Reshape to matrix: rows = regions, cols = cell types
  mat_df <- tidyr::pivot_wider(prop_df,
                               names_from = celltype_col,
                               values_from = proportion_col,
                               values_fill = 0)
  
  mat <- as.data.frame(mat_df)
  rownames(mat) <- mat[[region_col]]
  mat[[region_col]] <- NULL
  mat <- as.matrix(mat)
  
  # Scale matrix if requested
  if (identical(scale, "row")) {
    mat <- t(scale(t(mat)))
  } else if (identical(scale, "column")) {
    mat <- scale(mat)
  }
  
  # Plot heatmap
  pheatmap::pheatmap(mat,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     color = color_palette(100),
                     main = "Cell Type Proportions Heatmap")
  
  invisible(NULL)
}