#' Fit spatial GAM for gene expression across spatial coordinates
#'
#' Fits a GAM to model gene expression as a smooth function of spatial (x, y) coordinates.
#' Can fit models per group if `split_by` is provided.
#'
#' @param seurat_obj A Seurat object with spatial coordinates in metadata (e.g., "imagecol" and "imagerow").
#' @param gene Character; gene name to model.
#' @param coord_cols Character vector of length 2; names of metadata columns with spatial coords (default c("imagecol", "imagerow")).
#' @param assay Assay name to get expression data from (default "RNA").
#' @param slot Slot name of expression data in assay (default "data").
#' @param split_by Optional metadata column name to fit separate GAMs per group.
#' @param family Family object for GAM (default Gaussian).
#' @param k Number of basis functions for smooth term (default 10).
#'
#' @return A list containing:
#'   - gam_fits: list of fitted mgcv::gam models (one per group or single).
#'   - predictions: data.frame with predicted expression on a spatial grid per group.
#'   - grid: data.frame with spatial coordinates used for prediction.
#' @importFrom mgcv gam s
#' @export
#'
#' @examples
#' fit_res <- fit_spatial_gam(seurat_obj, gene = "ACTA2", split_by = "sample_id")
fit_spatial_gam <- function(seurat_obj,
                            gene,
                            coord_cols = c("imagecol", "imagerow"),
                            assay = "RNA",
                            slot = "data",
                            split_by = NULL,
                            family = gaussian(),
                            k = 10) {
  library(mgcv)
  meta <- seurat_obj@meta.data
  
  # Validate inputs
  if (!all(coord_cols %in% colnames(meta))) {
    stop(paste("Coordinate columns", paste(coord_cols, collapse = ", "), "not found in metadata."))
  }
  if (!gene %in% rownames(seurat_obj)) {
    stop(paste("Gene", gene, "not found in Seurat object."))
  }
  if (!is.null(split_by) && !split_by %in% colnames(meta)) {
    stop(paste("Split_by column", split_by, "not found in metadata."))
  }
  
  expr_vec <- as.numeric(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, ])
  df <- meta[, coord_cols, drop = FALSE]
  df$expression <- expr_vec
  if (!is.null(split_by)) df$group <- as.factor(meta[[split_by]]) else df$group <- factor("all")
  
  gam_fits <- list()
  predictions <- list()
  grid_list <- list()
  
  for (grp in levels(df$group)) {
    df_grp <- df[df$group == grp, , drop = FALSE]
    
    # Fit GAM: expression ~ s(x, y)
    gam_fit <- mgcv::gam(expression ~ s(get(coord_cols[1]), get(coord_cols[2]), k = k),
                         data = df_grp,
                         family = family)
    
    gam_fits[[grp]] <- gam_fit
    
    # Prepare prediction grid over the observed spatial range
    x_seq <- seq(min(df_grp[[coord_cols[1]]], na.rm = TRUE),
                 max(df_grp[[coord_cols[1]]], na.rm = TRUE), length.out = 100)
    y_seq <- seq(min(df_grp[[coord_cols[2]]], na.rm = TRUE),
                 max(df_grp[[coord_cols[2]]], na.rm = TRUE), length.out = 100)
    grid <- expand.grid(x = x_seq, y = y_seq)
    colnames(grid) <- coord_cols
    
    # Predict on grid
    grid$predicted_expression <- predict(gam_fit, newdata = grid, type = "response")
    grid$group <- grp
    
    predictions[[grp]] <- grid
    grid_list[[grp]] <- grid[, coord_cols, drop = FALSE]
  }
  
  predictions_df <- do.call(rbind, predictions)
  
  return(list(
    gam_fits = gam_fits,
    predictions = predictions_df,
    grids = grid_list
  ))
}