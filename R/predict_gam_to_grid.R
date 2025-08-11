#' Predict GAM smoothed expression values on a spatial grid
#'
#' Given a fitted GAM model, predicts values on a regular grid over spatial coordinates.
#'
#' @param gam_model A fitted GAM model object (e.g., from mgcv::gam).
#' @param x_range Numeric vector of length 2 giving min and max x-coordinates to cover.
#' @param y_range Numeric vector of length 2 giving min and max y-coordinates to cover.
#' @param grid_length Integer number of grid points per axis (default 100).
#' @param group_var Optional character; name of grouping variable to expand grid by (default NULL).
#' @param group_levels Optional vector of group levels to use if group_var is provided.
#'
#' @return A data.frame with columns:
#'   - x: grid x coordinate
#'   - y: grid y coordinate
#'   - predicted: GAM predicted value at (x,y)
#'   - group: grouping variable level (if applicable)
#' @export
#'
#' @examples
#' grid_pred <- predict_gam_to_grid(gam_model, x_range = c(0, 1000), y_range = c(0, 1000))
predict_gam_to_grid <- function(gam_model,
                                x_range,
                                y_range,
                                grid_length = 100,
                                group_var = NULL,
                                group_levels = NULL) {
  if (length(x_range) != 2 || length(y_range) != 2) {
    stop("x_range and y_range must be numeric vectors of length 2.")
  }
  if (!inherits(gam_model, "gam")) {
    stop("gam_model must be a fitted mgcv::gam object.")
  }
  
  # Create base grid
  x_seq <- seq(x_range[1], x_range[2], length.out = grid_length)
  y_seq <- seq(y_range[1], y_range[2], length.out = grid_length)
  
  # If group_var is NULL, just predict on 2D grid
  if (is.null(group_var)) {
    grid_df <- expand.grid(x = x_seq, y = y_seq)
    colnames(grid_df) <- names(gam_model$model)[1:2]  # Assume first 2 model variables are x,y coords
    
    # Predict on grid
    grid_df$predicted <- predict(gam_model, newdata = grid_df, type = "response")
  } else {
    # group_var is specified - create grid expanded by group levels
    if (is.null(group_levels)) {
      stop("group_levels must be provided when group_var is specified.")
    }
    # Build grid with group factor column
    grids <- lapply(group_levels, function(g) {
      df <- expand.grid(x = x_seq, y = y_seq)
      colnames(df) <- names(gam_model$model)[1:2]
      df[[group_var]] <- g
      return(df)
    })
    
    grid_df <- do.call(rbind, grids)
    
    # Ensure group_var is factor with correct levels if model expects that
    if (is.factor(gam_model$model[[group_var]])) {
      grid_df[[group_var]] <- factor(grid_df[[group_var]], levels = levels(gam_model$model[[group_var]]))
    }
    
    # Predict on full grid
    grid_df$predicted <- predict(gam_model, newdata = grid_df, type = "response")
  }
  
  return(grid_df)
}