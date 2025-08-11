#' Calculate spatial GAM gradients for genes or signatures with visualization
#'
#' Fits GAMs to model spatial expression variation for specified genes or gene signatures,
#' optionally split by grouping variable, and produces smoothed spatial predictions and contour plots.
#'
#' @param seurat_obj A Seurat object containing spatial coordinates and expression data.
#' @param genes Character vector of gene names to model (or a single gene).
#' @param signature_list Named list of gene sets (optional); if provided, computes average signature expression.
#' @param assay Assay name for expression data (default "RNA").
#' @param slot Slot in assay for expression data (default "data").
#' @param spatial_coords Character vector of length 2 specifying spatial coordinate columns in metadata (default c("imagecol", "imagerow")).
#' @param split_by Optional metadata column name to split data by groups for separate GAM fits.
#' @param grid_length Integer; number of points per spatial axis for prediction grid (default 100).
#' @param plot_results Logical; whether to produce GAM contour plots (default TRUE).
#' @param gam_formula Optional; custom GAM formula (default: expression ~ s(x, y, bs = "tp")).
#' @param gam_family Family object for GAM fit (default gaussian()).
#' @param contour_bins Integer; number of contour levels in plots (default 10).
#' @param contour_alpha Numeric alpha for contour fill (default 0.8).
#' @param fill_palette Color palette function for plots (default viridis).
#'
#' @return A list with elements:
#'   - seurat_obj: the input Seurat object (unchanged).
#'   - gam_fits: named list of fitted GAM models per gene and group.
#'   - grid_predictions: data.frame of GAM predicted values on grid per gene and group.
#'   - plots: named list of ggplot2 contour plots (if plot_results=TRUE).
#'   - metadata: data.frame summarizing statistics per gene and group.
#' @import mgcv ggplot2 viridis dplyr tidyr
#' @export
#'
#' @examples
#' res <- spatial_gam_gradient(seurat_obj, genes = c("ACTA2", "COL1A1"), split_by = "tumor_core_location")
spatial_gam_gradient <- function(seurat_obj,
                                 genes,
                                 signature_list = NULL,
                                 assay = "RNA",
                                 slot = "data",
                                 spatial_coords = c("imagecol", "imagerow"),
                                 split_by = NULL,
                                 grid_length = 100,
                                 plot_results = TRUE,
                                 gam_formula = NULL,
                                 gam_family = mgcv::gaussian(),
                                 contour_bins = 10,
                                 contour_alpha = 0.8,
                                 fill_palette = viridis::viridis) {
  library(mgcv)
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(tidyr)
  
  meta <- seurat_obj@meta.data
  if (!all(spatial_coords %in% colnames(meta))) {
    stop("Spatial coordinate columns not found in metadata.")
  }
  
  # Prepare expression matrix or signature scores
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Helper: compute average signature expression per cell if signature_list provided
  compute_signature <- function(gene_vec) {
    genes_present <- intersect(gene_vec, rownames(expr_mat))
    if (length(genes_present) == 0) {
      warning("No genes found for signature: ", paste(gene_vec, collapse = ", "))
      return(rep(NA_real_, ncol(expr_mat)))
    }
    Matrix::colMeans(expr_mat[genes_present, , drop = FALSE])
  }
  
  # Prepare data.frame with spatial coords and optionally split_by group
  df_meta <- meta[, spatial_coords, drop = FALSE]
  if (!is.null(split_by)) {
    if (!(split_by %in% colnames(meta))) {
      stop("split_by column not found in metadata.")
    }
    df_meta[[split_by]] <- meta[[split_by]]
  } else {
    df_meta[[split_by]] <- factor("all")
  }
  
  # Prepare list of expression vectors to model: genes + signatures
  expr_list <- list()
  
  # Add genes
  for (g in genes) {
    if (!(g %in% rownames(expr_mat))) {
      warning("Gene not found in assay: ", g)
      next
    }
    expr_list[[g]] <- expr_mat[g, ]
  }
  
  # Add signatures if provided
  if (!is.null(signature_list)) {
    for (sig_name in names(signature_list)) {
      sig_scores <- compute_signature(signature_list[[sig_name]])
      expr_list[[sig_name]] <- sig_scores
    }
  }
  
  # Remove any NULL or NA expression sets
  expr_list <- expr_list[!sapply(expr_list, function(x) all(is.na(x)))]
  
  # For results storage
  gam_fits <- list()
  grid_predictions <- list()
  plots <- list()
  metadata <- data.frame()
  
  # Function to fit GAM and predict on grid for one gene and one group
  fit_predict_gam <- function(df_sub, expr_values, x_col, y_col, gam_formula, gam_family, grid_length) {
    # Prepare modeling data
    mod_df <- df_sub
    mod_df$expr <- expr_values
    
    # Remove NAs
    valid <- !is.na(mod_df$expr) & complete.cases(mod_df[, c(x_col, y_col)])
    mod_df <- mod_df[valid, ]
    
    if (nrow(mod_df) < 20) {
      warning("Not enough data points to fit GAM.")
      return(NULL)
    }
    
    # Default formula
    if (is.null(gam_formula)) {
      formula_use <- as.formula("expr ~ s(get(x_col), get(y_col), bs = 'tp')")
    } else {
      formula_use <- gam_formula
    }
    
    # Fit GAM
    gam_fit <- tryCatch({
      mgcv::gam(formula_use, data = mod_df, family = gam_family)
    }, error = function(e) {
      warning("GAM fit failed: ", e$message)
      return(NULL)
    })
    
    if (is.null(gam_fit)) return(NULL)
    
    # Create prediction grid
    x_seq <- seq(min(mod_df[[x_col]], na.rm = TRUE),
                 max(mod_df[[x_col]], na.rm = TRUE),
                 length.out = grid_length)
    y_seq <- seq(min(mod_df[[y_col]], na.rm = TRUE),
                 max(mod_df[[y_col]], na.rm = TRUE),
                 length.out = grid_length)
    grid_df <- expand.grid(x = x_seq, y = y_seq)
    colnames(grid_df) <- c(x_col, y_col)
    
    # Predict on grid
    grid_df$predicted <- predict(gam_fit, newdata = grid_df, type = "response")
    
    return(list(gam_fit = gam_fit, grid_df = grid_df))
  }
  
  # Iterate over groups and genes/signatures
  groups <- unique(df_meta[[split_by]])
  for (grp in groups) {
    df_sub <- df_meta[df_meta[[split_by]] == grp, , drop = FALSE]
    
    for (nm in names(expr_list)) {
      expr_values <- expr_list[[nm]][colnames(seurat_obj) %in% rownames(df_sub)]
      
      # Align cell names
      cells_in_group <- intersect(colnames(seurat_obj), rownames(df_sub))
      expr_values <- expr_list[[nm]][cells_in_group]
      
      if (length(expr_values) == 0) next
      
      res <- fit_predict_gam(df_sub[cells_in_group, , drop = FALSE], expr_values,
                             spatial_coords[1], spatial_coords[2],
                             gam_formula, gam_family, grid_length)
      
      if (is.null(res)) next
      
      gam_fits[[paste(nm, grp, sep = "_")]] <- res$gam_fit
      grid_df <- res$grid_df
      grid_df[[split_by]] <- grp
      grid_df$feature <- nm
      grid_predictions[[paste(nm, grp, sep = "_")]] <- grid_df
      
      # Summary metadata for this fit
      r2 <- summary(res$gam_fit)$r.sq
      metadata <- rbind(metadata,
                        data.frame(feature = nm,
                                   group = grp,
                                   r_squared = r2,
                                   stringsAsFactors = FALSE))
    }
  }
  
  # Combine grid predictions into single data.frame
  grid_pred_df <- do.call(rbind, grid_predictions)
  
  # Plotting
  if (plot_results) {
    for (nm_grp in names(grid_predictions)) {
      p <- ggplot(grid_predictions[[nm_grp]], aes_string(x = spatial_coords[1], y = spatial_coords[2], fill = "predicted")) +
        geom_raster(interpolate = TRUE) +
        scale_fill_viridis(option = "magma", direction = -1, alpha = contour_alpha) +
        coord_fixed() +
        geom_contour(aes(z = predicted), color = "white", bins = contour_bins, size = 0.3) +
        theme_minimal() +
        labs(fill = "GAM predicted\nexpression",
             title = paste("GAM gradient:", nm_grp)) +
        theme(plot.title = element_text(hjust = 0.5))
      
      plots[[nm_grp]] <- p
    }
  }
  
  return(list(
    seurat_obj = seurat_obj,
    gam_fits = gam_fits,
    grid_predictions = grid_pred_df,
    plots = plots,
    metadata = metadata
  ))
}