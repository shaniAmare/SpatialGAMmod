#' Map contour region assignments to cells in a Seurat object
#'
#' Adds or updates a metadata column in the Seurat object that assigns each cell
#' to a spatial contour region based on a provided mapping data.frame.
#'
#' @param seurat_obj A Seurat object.
#' @param contour_mapping_df A data.frame with two columns:
#'    - cell: cell barcode or ID matching Seurat cell names.
#'    - contour: assigned contour/region label for each cell.
#' @param metadata_col Name of the metadata column to create/update (default "assigned_region").
#' @param overwrite Logical; whether to overwrite existing metadata column if it exists (default TRUE).
#'
#' @return The updated Seurat object with new/updated metadata column.
#' @export
#'
#' @examples
#' # contour_mapping_df with columns: cell, contour
#' seurat_obj <- map_contour_to_cells(seurat_obj, contour_mapping_df, metadata_col = "assigned_region")
map_contour_to_cells <- function(seurat_obj,
                                 contour_mapping_df,
                                 metadata_col = "assigned_region",
                                 overwrite = TRUE) {
  meta <- seurat_obj@meta.data
  
  # Validate inputs
  if (!all(c("cell", "contour") %in% colnames(contour_mapping_df))) {
    stop("contour_mapping_df must contain columns named 'cell' and 'contour'.")
  }
  
  cells_in_seurat <- rownames(meta)
  mapping_cells <- contour_mapping_df$cell
  
  # Filter mapping to cells present in Seurat object
  valid_mapping <- contour_mapping_df[mapping_cells %in% cells_in_seurat, ]
  
  if (nrow(valid_mapping) == 0) {
    stop("No matching cells found between contour_mapping_df and Seurat object.")
  }
  
  # Prepare vector for assignment
  new_col <- rep(NA_character_, length(cells_in_seurat))
  names(new_col) <- cells_in_seurat
  
  # Assign contour labels from mapping
  new_col[valid_mapping$cell] <- as.character(valid_mapping$contour)
  
  # Check if metadata_col exists
  if (metadata_col %in% colnames(meta) && !overwrite) {
    warning(paste0("Metadata column '", metadata_col, "' exists and overwrite=FALSE. No changes made."))
    return(seurat_obj)
  }
  
  # Add or update metadata column
  meta[[metadata_col]] <- new_col
  
  seurat_obj@meta.data <- meta
  
  return(seurat_obj)
}