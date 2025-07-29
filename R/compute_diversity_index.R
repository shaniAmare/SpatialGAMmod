#' Compute cell type diversity index per spatial region
#'
#' This function calculates ecological diversity metrics (e.g., Shannon, Simpson)
#' for cell type distributions within spatial regions (e.g., GAM contours).
#'
#' @param seurat_obj A Seurat object with metadata columns for region and cell type.
#' @param region_col Metadata column indicating spatial region or contour (e.g., "assigned_region").
#' @param celltype_col Metadata column indicating cell type annotations.
#' @param index Diversity index to compute. One of: "shannon" (default), "simpson", or "invsimpson".
#'
#' @return A data.frame with columns: region, diversity_index.
#' @export
#'
#' @examples
#' diversity_df <- compute_diversity_index(seurat_obj,
#'                                         region_col = "assigned_region",
#'                                         celltype_col = "CellType",
#'                                         index = "shannon")
compute_diversity_index <- function(seurat_obj,
                                    region_col = "assigned_region",
                                    celltype_col = "CellType",
                                    index = c("shannon", "simpson", "invsimpson")) {
  index <- match.arg(index)
  
  meta <- seurat_obj@meta.data
  
  if (!region_col %in% colnames(meta)) {
    stop(paste("Region column", region_col, "not found in metadata."))
  }
  
  if (!celltype_col %in% colnames(meta)) {
    stop(paste("Cell type column", celltype_col, "not found in metadata."))
  }
  
  # Count cell types per region
  count_df <- meta |>
    dplyr::count(.data[[region_col]], .data[[celltype_col]]) |>
    dplyr::rename(region = 1, celltype = 2, count = 3)
  
  # Convert to wide matrix (rows = regions, columns = celltypes)
  count_matrix <- tidyr::pivot_wider(count_df,
                                     names_from = celltype,
                                     values_from = count,
                                     values_fill = 0) |>
    tibble::column_to_rownames("region") |>
    as.matrix()
  
  # Compute diversity
  diversity_vals <- vegan::diversity(count_matrix, index = index)
  
  result_df <- data.frame(
    region = names(diversity_vals),
    diversity_index = as.numeric(diversity_vals),
    row.names = NULL
  )
  
  return(result_df)
}
