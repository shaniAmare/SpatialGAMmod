#' Assess variability of values across contour regions grouped by a condition
#'
#' This function summarizes variability of a feature (e.g., gene expression or cell type proportion)
#' within each contour region, comparing values across different groupings.
#'
#' @param df A data frame containing at least these columns: region, group, value.
#' @param region_col Name of the column indicating spatial region (e.g., "assigned_region").
#' @param group_col Name of the column indicating group/condition (e.g., "TissueType").
#' @param value_col Name of the column with numeric values to summarize.
#' @param stat_function Statistical function to compute variability (default: `sd` for standard deviation).
#'
#' @return A data.frame with region-wise variability statistics across groups.
#' @export
#'
#' @examples
#' variability_df <- contour_variability_by_group(
#'   df = prop_df,
#'   region_col = "assigned_region",
#'   group_col = "TissueType",
#'   value_col = "proportion",
#'   stat_function = sd
#' )
contour_variability_by_group <- function(df,
                                         region_col = "assigned_region",
                                         group_col = "group",
                                         value_col = "value",
                                         stat_function = sd) {
  # Check columns
  required_cols <- c(region_col, group_col, value_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing one or more required columns in input data frame.")
  }
  
  df_summary <- df |>
    dplyr::group_by(.data[[region_col]]) |>
    dplyr::summarise(
      n_groups = dplyr::n_distinct(.data[[group_col]]),
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      variability = stat_function(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_summary)
}