#' Plot spatial GAM smoothed expression as contour maps
#'
#' Visualizes GAM predicted spatial expression values on spatial coordinates as contour plots.
#'
#' @param gam_pred_df A data.frame containing columns:
#'   - x: spatial x coordinate (numeric)
#'   - y: spatial y coordinate (numeric)
#'   - predicted: GAM predicted expression or signature value (numeric)
#'   - group: optional grouping variable (factor or character) for faceting or coloring (optional)
#' @param color_palette Color palette function or vector for fill colors (default viridis).
#' @param facet_by Character; column name in `gam_pred_df` to facet plots by (default NULL, no faceting).
#' @param point_size Numeric size of points underlying contours (default 0.5).
#' @param contour_lines Logical; whether to add contour lines (default TRUE).
#' @param contour_bins Integer; number of contour bins/levels (default 10).
#' @param alpha_fill Numeric transparency for fill colors (default 0.8).
#'
#' @return A ggplot2 object of the contour plot(s).
#' @import ggplot2
#' @export
#'
#' @examples
#' plot_gam_contours(gam_pred_df, facet_by = "group")
plot_gam_contours <- function(gam_pred_df,
                              color_palette = viridis::viridis,
                              facet_by = NULL,
                              point_size = 0.5,
                              contour_lines = TRUE,
                              contour_bins = 10,
                              alpha_fill = 0.8) {
  library(ggplot2)
  library(viridis)
  
  required_cols <- c("x", "y", "predicted")
  if (!all(required_cols %in% colnames(gam_pred_df))) {
    stop("gam_pred_df must contain columns: x, y, predicted")
  }
  
  p <- ggplot(gam_pred_df, aes(x = x, y = y))
  
  if (!is.null(facet_by) && facet_by %in% colnames(gam_pred_df)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  
  # Use fill for predicted values
  p <- p +
    geom_raster(aes(fill = predicted), interpolate = TRUE) +
    scale_fill_viridis(option = "magma", direction = -1, alpha = alpha_fill) +
    coord_fixed() +
    labs(fill = "GAM predicted\nexpression") +
    theme_minimal()
  
  if (contour_lines) {
    p <- p + geom_contour(aes(z = predicted), bins = contour_bins, color = "white", size = 0.3)
  }
  
  # Optionally add points of underlying grid with low alpha
  if ("group" %in% colnames(gam_pred_df)) {
    p <- p + geom_point(aes(color = group), size = point_size, alpha = 0.6) +
      guides(color = guide_legend(title = "Group"))
  } else {
    p <- p + geom_point(size = point_size, color = "black", alpha = 0.3)
  }
  
  return(p)
}