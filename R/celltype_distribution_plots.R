#' Plot cell type distribution across spatial regions
#'
#' This function visualizes the distribution of cell types across spatial regions
#' as either a barplot (stacked or filled) or a heatmap.
#'
#' @param proportions_df A data.frame from `calculate_celltype_proportions()`.
#'                       Must contain columns: region, celltype, proportion or count.
#' @param plot_type One of "bar" (default), "filled_bar", or "heatmap".
#' @param use_counts Logical; if TRUE, uses `count` instead of `proportion` (if available).
#' @param palette Optional color palette for cell types.
#' @param region_order Optional vector to control the order of regions on the x-axis or y-axis.
#' @param celltype_order Optional vector to control the order of cell types.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' proportions_df <- calculate_celltype_proportions(seurat_obj)
#' celltype_distribution_plots(proportions_df, plot_type = "bar")
celltype_distribution_plots <- function(proportions_df,
                                        plot_type = c("bar", "filled_bar", "heatmap"),
                                        use_counts = FALSE,
                                        palette = NULL,
                                        region_order = NULL,
                                        celltype_order = NULL) {
  plot_type <- match.arg(plot_type)
  
  df <- proportions_df
  
  # Check required columns
  required_cols <- c("region", "celltype", if (use_counts) "count" else "proportion")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  value_col <- if (use_counts) "count" else "proportion"
  
  # Factor ordering
  if (!is.null(region_order)) {
    df$region <- factor(df$region, levels = region_order)
  }
  
  if (!is.null(celltype_order)) {
    df$celltype <- factor(df$celltype, levels = celltype_order)
  }
  
  if (plot_type %in% c("bar", "filled_bar")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = region, y = .data[[value_col]], fill = celltype)) +
      ggplot2::geom_bar(stat = "identity", position = ifelse(plot_type == "filled_bar", "fill", "stack")) +
      ggplot2::labs(x = "Region", y = ifelse(use_counts, "Cell Count", "Proportion"),
                    fill = "Cell Type") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (!is.null(palette)) {
      p <- p + ggplot2::scale_fill_manual(values = palette)
    }
    
  } else if (plot_type == "heatmap") {
    # Fill in missing combinations with zero
    df_complete <- tidyr::complete(df, region, celltype, fill = list(!!value_col := 0))
    
    p <- ggplot2::ggplot(df_complete, ggplot2::aes(x = region, y = celltype, fill = .data[[value_col]])) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
      ggplot2::labs(x = "Region", y = "Cell Type", fill = ifelse(use_counts, "Count", "Proportion")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}