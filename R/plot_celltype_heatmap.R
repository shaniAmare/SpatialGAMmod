#' Plot heatmap of cell type proportions or counts across regions
#'
#' Creates a heatmap showing cell type proportions or counts per spatial region/group.
#'
#' @param df Data frame with columns: region/group, cell type, and numeric value (proportion or count).
#' @param region_col Name of the column for regions/groups (default "assigned_region").
#' @param celltype_col Name of the column for cell type labels (default "CellType").
#' @param value_col Name of the numeric column for proportion/count (default "proportion").
#' @param scale Logical; whether to scale values per row or column ("row", "column", or FALSE) (default FALSE).
#' @param cluster_rows Logical; whether to cluster rows (default TRUE).
#' @param cluster_cols Logical; whether to cluster columns (default TRUE).
#' @param color_palette Color palette function or vector (default viridis).
#'
#' @return A ggplot2 heatmap object or a `pheatmap` plot (depending on availability).
#' @import ggplot2 tidyr dplyr
#' @export
#'
#' @examples
#' # df must have columns assigned_region, CellType, proportion
#' plot_celltype_heatmap(df, region_col = "assigned_region", celltype_col = "CellType")
plot_celltype_heatmap <- function(df,
                                  region_col = "assigned_region",
                                  celltype_col = "CellType",
                                  value_col = "proportion",
                                  scale = FALSE,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  color_palette = viridis::viridis) {
  # Try to load pheatmap if available; else fallback to ggplot2 heatmap
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    # Prepare matrix: rows=regions, cols=celltypes
    mat_df <- df %>%
      tidyr::pivot_wider(names_from = celltype_col, values_from = value_col, values_fill = 0) %>%
      tibble::column_to_rownames(region_col)
    
    mat <- as.matrix(mat_df)
    
    # Scale matrix if requested
    if (identical(scale, "row")) {
      mat <- t(scale(t(mat)))
    } else if (identical(scale, "column")) {
      mat <- scale(mat)
    }
    
    pheatmap::pheatmap(mat,
                       cluster_rows = cluster_rows,
                       cluster_cols = cluster_cols,
                       color = color_palette(100))
  } else {
    # Use ggplot2 heatmap
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    
    # Optionally scale within df before plotting
    df_plot <- df
    if (identical(scale, "row")) {
      df_plot <- df_plot %>%
        group_by(!!rlang::sym(region_col)) %>%
        mutate(scaled_value = scale(!!rlang::sym(value_col))) %>%
        ungroup()
      val_col <- "scaled_value"
    } else if (identical(scale, "column")) {
      df_plot <- df_plot %>%
        group_by(!!rlang::sym(celltype_col)) %>%
        mutate(scaled_value = scale(!!rlang::sym(value_col))) %>%
        ungroup()
      val_col <- "scaled_value"
    } else {
      val_col <- value_col
    }
    
    p <- ggplot(df_plot, aes_string(x = celltype_col, y = region_col, fill = val_col)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab(celltype_col) +
      ylab(region_col) +
      labs(fill = value_col)
    
    print(p)
  }
}