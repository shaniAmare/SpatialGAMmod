#' Plot expression density distributions per spatial region
#'
#' Visualizes the distribution of expression values for a gene or signature across contour regions.
#'
#' @param seurat_obj A Seurat object with expression and metadata.
#' @param gene Character; gene name to plot.
#' @param contour_col Metadata column specifying contour/region assignment (default "assigned_region").
#' @param assay Assay name to use for expression (default "RNA").
#' @param slot Slot to use in assay (default "data").
#' @param plot_type Character; one of "violin", "boxplot", or "density" (default "violin").
#' @param log_scale Logical; if TRUE, log-transform expression (log1p) before plotting (default TRUE).
#' @param fill_palette Optional vector or palette function for fill colors by region.
#'
#' @return A ggplot2 plot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' plot_expression_density_per_region(seurat_obj, gene = "ACTA2", contour_col = "assigned_region")
plot_expression_density_per_region <- function(seurat_obj,
                                               gene,
                                               contour_col = "assigned_region",
                                               assay = "RNA",
                                               slot = "data",
                                               plot_type = c("violin", "boxplot", "density"),
                                               log_scale = TRUE,
                                               fill_palette = NULL) {
  library(ggplot2)
  
  plot_type <- match.arg(plot_type)
  meta <- seurat_obj@meta.data
  
  if (!(contour_col %in% colnames(meta))) {
    stop(paste0("Metadata column '", contour_col, "' not found in Seurat object."))
  }
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste0("Gene '", gene, "' not found in Seurat object."))
  }
  
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, ]
  df <- data.frame(
    expression = as.numeric(expr),
    region = factor(meta[[contour_col]], levels = unique(meta[[contour_col]]))
  )
  
  if (log_scale) {
    df$expression <- log1p(df$expression)
  }
  
  p <- ggplot(df, aes(x = region, y = expression, fill = region))
  
  if (plot_type == "violin") {
    p <- p + geom_violin(trim = TRUE, scale = "width")
  } else if (plot_type == "boxplot") {
    p <- p + geom_boxplot()
  } else if (plot_type == "density") {
    # Density plot per region requires facetting or color by region
    p <- ggplot(df, aes(x = expression, color = region, fill = region)) +
      geom_density(alpha = 0.4) +
      theme_minimal() +
      labs(x = ifelse(log_scale, paste0("log1p(", gene, ")"), gene),
           y = "Density") +
      theme(legend.title = element_blank())
    return(p)
  }
  
  p <- p +
    theme_minimal() +
    xlab(contour_col) +
    ylab(ifelse(log_scale, paste0("log1p(", gene, ")"), gene)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = FALSE)
  
  if (!is.null(fill_palette)) {
    p <- p + scale_fill_manual(values = fill_palette)
  } else {
    p <- p + scale_fill_brewer(palette = "Set3")
  }
  
  return(p)
}