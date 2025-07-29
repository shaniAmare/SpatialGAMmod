#' Compare gene expression across spatial contour regions
#'
#' This function tests for differential expression of genes or gene signatures across
#' spatial regions (e.g., contour assignments), using either pairwise or group-wise comparisons.
#'
#' @param seurat_obj A Seurat object with spatial region annotations and expression data.
#' @param genes Character vector of gene names, or a named list of gene sets to average.
#' @param region_col Metadata column indicating spatial region or contour groupings.
#' @param assay Assay to pull expression data from. Defaults to "RNA".
#' @param slot Slot to extract expression values from. Defaults to "data".
#' @param test_method Statistical test to use: "wilcox" (default), "anova", or "kruskal".
#' @param min_cells Minimum number of cells per region to include in test.
#' @param p_adjust_method Multiple testing correction method for p-values. Default: "fdr".
#'
#' @return A data.frame with expression comparisons per gene/signature across regions.
#' @export
#'
#' @examples
#' compare_df <- compare_expression_across_contours(seurat_obj,
#'                    genes = c("ACTA2", "COL1A1"),
#'                    region_col = "assigned_region")
compare_expression_across_contours <- function(seurat_obj,
                                               genes,
                                               region_col = "assigned_region",
                                               assay = "RNA",
                                               slot = "data",
                                               test_method = c("wilcox", "anova", "kruskal"),
                                               min_cells = 10,
                                               p_adjust_method = "fdr") {
  test_method <- match.arg(test_method)
  meta <- seurat_obj@meta.data
  
  if (!region_col %in% colnames(meta)) {
    stop(paste("Region column", region_col, "not found in metadata."))
  }
  
  expr_data <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  region_vector <- meta[[region_col]]
  
  # Prepare result container
  results_list <- list()
  
  # If genes is a list (e.g., gene signatures), compute mean signature scores
  if (is.list(genes)) {
    gene_scores <- sapply(names(genes), function(sig_name) {
      gene_set <- genes[[sig_name]]
      valid_genes <- intersect(gene_set, rownames(expr_data))
      if (length(valid_genes) == 0) return(rep(NA, ncol(seurat_obj)))
      Matrix::colMeans(expr_data[valid_genes, , drop = FALSE])
    })
    gene_scores <- as.data.frame(gene_scores)
    colnames(gene_scores) <- names(genes)
  } else {
    valid_genes <- intersect(genes, rownames(expr_data))
    if (length(valid_genes) == 0) stop("None of the requested genes were found.")
    gene_scores <- as.data.frame(t(as.matrix(expr_data[valid_genes, , drop = FALSE])))
  }
  
  gene_scores$region <- region_vector
  
  for (gene in setdiff(colnames(gene_scores), "region")) {
    df <- gene_scores[, c("region", gene), drop = FALSE]
    colnames(df) <- c("region", "expr")
    df <- df[!is.na(df$region) & !is.na(df$expr), ]
    region_counts <- table(df$region)
    df <- df[df$region %in% names(region_counts[region_counts >= min_cells]), ]
    
    if (nrow(df) == 0 || length(unique(df$region)) < 2) {
      next
    }
    
    # Perform the chosen test
    test_result <- tryCatch({
      if (test_method == "wilcox") {
        stats::pairwise.wilcox.test(df$expr, df$region, p.adjust.method = p_adjust_method)
      } else if (test_method == "anova") {
        aov_out <- stats::aov(expr ~ region, data = df)
        summary(aov_out)[[1]][["Pr(>F)"]][1]
      } else if (test_method == "kruskal") {
        stats::kruskal.test(expr ~ region, data = df)$p.value
      }
    }, error = function(e) NA)
    
    if (test_method == "wilcox") {
      pval_mat <- test_result$p.value
      pval_df <- as.data.frame(as.table(pval_mat))
      colnames(pval_df) <- c("group1", "group2", "p_value")
      pval_df$gene <- gene
      results_list[[length(results_list) + 1]] <- pval_df
    } else {
      results_list[[length(results_list) + 1]] <- data.frame(
        gene = gene,
        test = test_method,
        p_value = test_result,
        stringsAsFactors = FALSE
      )
    }
  }
  
  results_df <- do.call(rbind, results_list)
  
  # Handle p-adjust if not done already (e.g., for ANOVA or Kruskal)
  if (!"p_adj" %in% colnames(results_df)) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  }
  
  return(results_df)
}
