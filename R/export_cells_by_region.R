#' Export cell metadata and optionally expression data by spatial region
#'
#' Extracts cells from specified spatial regions and exports their metadata and/or
#' expression data to CSV files or returns as a named list of data.frames.
#'
#' @param seurat_obj A Seurat object containing spatial region metadata.
#' @param region_col Metadata column indicating spatial region assignments (e.g., "assigned_region").
#' @param regions Character vector specifying which regions to export. Default is all regions.
#' @param export_expression Logical; whether to include expression matrix for exported cells. Default FALSE.
#' @param assay Assay name to extract expression from (default "RNA").
#' @param slot Slot of assay data to extract (default "counts").
#' @param output_dir Optional directory path to save CSV files. If NULL, no files are saved.
#' @param prefix Optional filename prefix for saved files.
#'
#' @return Named list of lists: each element named by region, containing:
#'   - metadata: data.frame of cell metadata
#'   - expression: sparse or dense matrix of expression (if export_expression = TRUE)
#' @export
#'
#' @examples
#' # Export all regions metadata only
#' export_cells_by_region(seurat_obj, region_col = "assigned_region", output_dir = "exports/")
export_cells_by_region <- function(seurat_obj,
                                   region_col = "assigned_region",
                                   regions = NULL,
                                   export_expression = FALSE,
                                   assay = "RNA",
                                   slot = "counts",
                                   output_dir = NULL,
                                   prefix = "cells") {
  meta <- seurat_obj@meta.data
  
  if (!region_col %in% colnames(meta)) {
    stop(paste("Region column", region_col, "not found in metadata."))
  }
  
  if (is.null(regions)) {
    regions <- unique(meta[[region_col]])
  } else {
    missing_regions <- setdiff(regions, unique(meta[[region_col]]))
    if (length(missing_regions) > 0) {
      warning("The following regions are not present in the data and will be skipped: ",
              paste(missing_regions, collapse = ", "))
      regions <- intersect(regions, unique(meta[[region_col]]))
    }
  }
  
  results <- list()
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  for (region in regions) {
    cells_in_region <- rownames(meta)[meta[[region_col]] == region]
    metadata_region <- meta[cells_in_region, , drop = FALSE]
    
    export_list <- list(metadata = metadata_region)
    
    if (export_expression) {
      expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
      expr_region <- expr_mat[, cells_in_region, drop = FALSE]
      export_list$expression <- expr_region
    }
    
    results[[region]] <- export_list
    
    # Write files if output_dir specified
    if (!is.null(output_dir)) {
      meta_file <- file.path(output_dir, paste0(prefix, "_", region, "_metadata.csv"))
      utils::write.csv(metadata_region, meta_file, row.names = TRUE)
      
      if (export_expression) {
        expr_file <- file.path(output_dir, paste0(prefix, "_", region, "_expression.csv"))
        # For expression, write dense matrix (may be large)
        expr_df <- as.data.frame(as.matrix(export_list$expression))
        utils::write.csv(expr_df, expr_file, row.names = TRUE)
      }
    }
  }
  
  return(results)
}