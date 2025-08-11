#' Compare gene expression between tumour and non-tumour contour regions
#'
#' Performs differential expression/enrichment analysis between tumour vs non-tumour regions
#' based on contour assignments in metadata, and plots the results.
#'
#' @param seurat_obj A Seurat object with spatial metadata containing a contour region annotation.
#' @param contour_col Metadata column defining contour regions (e.g., "assigned_region").
#' @param tumour_regions Character vector of contour region names considered tumour (e.g., c("tumour_core", "tumour_invasive")).
#' @param non_tumour_regions Character vector of contour region names considered non-tumour (e.g., c("stroma", "immune")).
#' @param genes Character vector of gene names to test (default NULL means all genes in assay).
#' @param assay Assay name for expression data (default "RNA").
#' @param slot Slot in assay to use (default "data").
#' @param test_method Statistical test to use (default "wilcox").
#' @param p_adjust_method Method for multiple testing correction (default "BH").
#' @param plot_results Logical; whether to generate a volcano plot (default TRUE).
#'
#' @return A list containing:
#'   - results: data.frame with differential expression statistics per gene.
#'   - plot: volcano plot (ggplot) if plot_results = TRUE, else NULL.
#' @importFrom Seurat GetAssayData
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline scale_color_manual theme_minimal labs
#' @export
#'
#' @examples
#' res <- tumour_vs_non_tumour_contour(seurat_obj,
#'                                    contour_col = "assigned_region",
#'                                    tumour_regions = c("tumour_core"),
#'                                    non_tumour_regions = c("stroma", "immune"),
#'                                    genes = c("ACTA2", "COL1A1"))
tumour_vs_non_tumour_contour <- function(seurat_obj,
                                         contour_col,
                                         tumour_regions,
                                         non_tumour_regions,
                                         genes = NULL,
                                         assay = "RNA",
                                         slot = "data",
                                         test_method = "wilcox",
                                         p_adjust_method = "BH",
                                         plot_results = TRUE) {
  if (!(contour_col %in% colnames(seurat_obj@meta.data))) {
    stop(paste0("Contour column '", contour_col, "' not found in Seurat metadata."))
  }
  
  meta <- seurat_obj@meta.data
  cells_tumour <- rownames(meta)[meta[[contour_col]] %in% tumour_regions]
  cells_non_tumour <- rownames(meta)[meta[[contour_col]] %in% non_tumour_regions]
  
  if (length(cells_tumour) < 3 || length(cells_non_tumour) < 3) {
    stop("Not enough cells in tumour or non-tumour groups for testing.")
  }
  
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Use all genes if none specified
  if (is.null(genes)) {
    genes <- rownames(expr_mat)
  } else {
    genes <- intersect(genes, rownames(expr_mat))
  }
  
  # Prepare results storage
  results <- data.frame(
    gene = genes,
    p_value = NA_real_,
    avg_log2FC = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Compute avg log2 fold change and wilcox p-values
  for (i in seq_along(genes)) {
    g <- genes[i]
    
    expr_tumour <- expr_mat[g, cells_tumour]
    expr_non_tumour <- expr_mat[g, cells_non_tumour]
    
    # Compute average log2 fold change
    avg_expr_tumour <- mean(expr_tumour, na.rm = TRUE)
    avg_expr_non_tumour <- mean(expr_non_tumour, na.rm = TRUE)
    avg_log2FC <- log2(avg_expr_tumour + 1e-6) - log2(avg_expr_non_tumour + 1e-6)
    results$avg_log2FC[i] <- avg_log2FC
    
    # Statistical test
    if (test_method == "wilcox") {
      test_res <- wilcox.test(expr_tumour, expr_non_tumour)
      p_val <- test_res$p.value
    } else {
      stop("Only 'wilcox' test_method currently supported.")
    }
    results$p_value[i] <- p_val
  }
  
  # Adjust p-values
  results$p_adj <- p.adjust(results$p_value, method = p_adjust_method)
  
  # Order by significance
  results <- results[order(results$p_adj), ]
  
  # Volcano plot
  plot_obj <- NULL
  if (plot_results) {
    library(ggplot2)
    results$significant <- results$p_adj < 0.05
    results$log10p <- -log10(results$p_adj + 1e-10)
    
    plot_obj <- ggplot(results, aes(x = avg_log2FC, y = log10p, color = significant)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("grey", "red")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_minimal() +
      labs(title = "Tumour vs Non-Tumour Differential Expression",
           x = "Average Log2 Fold Change",
           y = "-log10 Adjusted p-value",
           color = "Significant")
  }
  
  return(list(results = results, plot = plot_obj))
}