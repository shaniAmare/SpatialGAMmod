#' Compute Shannon entropy of cell types per spatial region
#'
#' This function computes the Shannon entropy of cell type distributions
#' within each spatial region (e.g., GAM contour), based on metadata in a Seurat object.
#'
#' @param seurat_obj A Seurat object with region and cell type metadata.
#' @param region_col Name of the metadata column for region assignment.
#' @param celltype_col Name of the metadata column for cell type annotation.
#' @param base Base of logarithm used for entropy. Default is 2 (binary entropy).
#'
#' @return A data.frame with columns: region, entropy.
#' @export
#'
#' @examples
#' entropy_df <- compute_entropy_per_region(seurat_obj,
#'                                          region_col = "assigned_region",
#'                                          celltype_col = "CellType")
compute_entropy_per_region <- function(seurat_obj,
                                       region_col = "assigned_region",
                                       celltype_col = "CellType",
                                       base = 2) {
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
  
  # Normalize to proportions and compute entropy per region
  entropy_df <- count_df |>
    dplyr::group_by(region) |>
    dplyr::mutate(prop = count / sum(count)) |>
    dplyr::summarise(
      entropy = -sum(prop * log(prop, base = base), na.rm = TRUE),
      .groups = "drop"
    )
  
  return(entropy_df)
}