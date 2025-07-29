#' Assign spatial GAM contour regions to cells
#'
#' This function assigns GAM-derived contour regions (e.g., isocontours) to cells
#' in a Seurat object based on a precomputed metadata column (e.g., "contour_id").
#' If a `split_by` column is provided, region labels are made unique across groups.
#'
#' @param seurat_obj A Seurat object containing GAM-derived contour assignments.
#' @param contour_col Name of the metadata column containing the raw contour IDs (e.g., "gam_contour").
#' @param new_col Name of the new column to store the final assigned region labels. Defaults to "assigned_region".
#' @param split_by Optional metadata column to split by (e.g., "tumor_core_location" or "Run_ID").
#'        Region labels will be prefixed with the split group for uniqueness.
#' @param drop_na Logical; if TRUE, cells with NA contours are dropped. If FALSE, labeled as "Unknown".
#'
#' @return The input Seurat object with an additional column in `@meta.data` for assigned region.
#' @export
#'
#' @examples
#' seurat_obj <- assign_contour_regions(seurat_obj, contour_col = "gam_contour", split_by = "sample_id")
assign_contour_regions <- function(seurat_obj,
                                   contour_col = "gam_contour",
                                   new_col = "assigned_region",
                                   split_by = NULL,
                                   drop_na = FALSE) {
  meta <- seurat_obj@meta.data
  
  if (!contour_col %in% colnames(meta)) {
    stop(paste("Contour column", contour_col, "not found in metadata."))
  }
  
  if (!is.null(split_by) && !split_by %in% colnames(meta)) {
    stop(paste("Split column", split_by, "not found in metadata."))
  }
  
  # Base assignment
  region_vec <- meta[[contour_col]]
  
  # Handle split_by
  if (!is.null(split_by)) {
    group_vec <- as.character(meta[[split_by]])
    region_vec <- ifelse(is.na(region_vec), NA,
                         paste0(group_vec, "_region_", region_vec))
  } else {
    region_vec <- ifelse(is.na(region_vec), NA,
                         paste0("region_", region_vec))
  }
  
  # Handle NA values
  if (drop_na) {
    region_vec <- region_vec[!is.na(region_vec)]
    seurat_obj <- subset(seurat_obj, cells = names(region_vec))
  } else {
    region_vec[is.na(region_vec)] <- "Unknown"
  }
  
  seurat_obj[[new_col]] <- region_vec
  
  return(seurat_obj)
}
