#' Calculate cell type proportions per spatial region
#'
#' This function computes the relative proportions of each cell type within each spatial region
#' (e.g., contour, cluster, core), based on metadata annotations in a Seurat object.
#'
#' @param seurat_obj A Seurat object with metadata columns for region and cell type.
#' @param region_col Name of the metadata column containing region assignments (e.g., "assigned_region").
#' @param celltype_col Name of the metadata column containing cell type annotations.
#' @param normalize Logical; if TRUE, proportions are returned (default). If FALSE, returns raw counts.
#'
#' @return A long-format `data.frame` with columns: region, celltype, count/proportion.
#' @export
#'
#' @examples
#' df <- calculate_celltype_proportions(seurat_obj,
#'                                      region_col = "assigned_region",
#'                                      celltype_col = "CellType")
calculate_celltype_proportions <- function(seurat_obj,
                                           region_col = "assigned_region",
                                           celltype_col = "CellType",
                                           normalize = TRUE) {
  meta <- seurat_obj@meta.data
  
  if (!region_col %in% colnames(meta)) {
    stop(paste("Region column", region_col, "not found in metadata."))
  }
  
  if (!celltype_col %in% colnames(meta)) {
    stop(paste("Cell type column", celltype_col, "not found in metadata."))
  }
  
  # Tabulate counts
  count_df <- meta |>
    dplyr::count(.data[[region_col]], .data[[celltype_col]]) |>
    dplyr::rename(region = 1, celltype = 2, count = 3)
  
  if (normalize) {
    count_df <- count_df |>
      dplyr::group_by(region) |>
      dplyr::mutate(proportion = count / sum(count)) |>
      dplyr::ungroup()
  }
  
  return(count_df)
}