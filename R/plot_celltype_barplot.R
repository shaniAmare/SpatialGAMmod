#' Plot cell type proportions or counts as barplots across regions or groups
#'
#' Generates a stacked or grouped barplot of cell type proportions or counts per contour region or group.
#'
#' @param df Data frame with columns: region/group, cell type, and count or proportion.
#' @param region_col Name of the column for spatial region or group (e.g., "assigned_region").
#' @param celltype_col Name of the column for cell type labels (e.g., "CellType").
#' @param value_col Name of the numeric column for counts or proportions (e.g., "proportion" or "count").
#' @param proportion Logical; if TRUE, expects proportions, else counts (default TRUE).
#' @param stacked Logical; if TRUE, plots stacked bars; if FALSE, grouped bars (default TRUE).
#' @param fill_palette Optional named vector or palette function for cell type colors.
#' @param title Plot title (optional).
#'
#' @return A ggplot2 barplot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # df should have columns: assigned_region, CellType, proportion
#' p <- plot_celltype_barplot(df, region_col = "assigned_region", celltype_col = "CellType", value_col = "proportion")
plot_celltype_barplot <- function(df,
                                  region_col = "assigned_region",
                                  celltype_col = "CellType",
                                  value_col = "proportion",
                                  proportion = TRUE,
                                  stacked = TRUE,
                                  fill_palette = NULL,
                                  title = NULL) {
  library(ggplot2)
  
  # Check required columns
  req_cols <- c(region_col, celltype_col, value_col)
  if (!all(req_cols %in% colnames(df))) {
    stop("Data frame must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  df[[region_col]] <- as.factor(df[[region_col]])
  df[[celltype_col]] <- as.factor(df[[celltype_col]])
  
  p <- ggplot(df, aes_string(x = region_col, y = value_col, fill = celltype_col))
  
  if (stacked) {
    p <- p + geom_bar(stat = "identity", position = "stack")
  } else {
    p <- p + geom_bar(stat = "identity", position = "dodge")
  }
  
  p <- p +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab(ifelse(proportion, "Proportion", "Count")) +
    xlab(region_col) +
    guides(fill = guide_legend(title = celltype_col))
  
  if (!is.null(fill_palette)) {
    p <- p + scale_fill_manual(values = fill_palette)
  } else {
    p <- p + scale_fill_brewer(palette = "Set3")
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}