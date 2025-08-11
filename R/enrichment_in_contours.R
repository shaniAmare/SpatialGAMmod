#' Enrichment analysis of gene(s) or signature(s) in spatial contours
#'
#' This function tests whether expression of specified genes or gene signatures
#' is significantly enriched in particular spatial contour regions compared to all other regions.
#'
#' @param seurat_obj A Seurat object with spatial contour metadata and expression data.
#' @param genes Character vector of genes, or named list of gene sets for signature scores.
#' @param contour_col Metadata column indicating contour/region assignments (e.g., "assigned_region").
#' @param assay Assay to use for expression (default "RNA").
#' @param slot Slot of assay data to use (default "data").
#' @param test_method Statistical test to compare expression between in-contour and out-contour cells.
#'                    Options: "wilcox" (default), "t.test".
#' @param min_cells Minimum number of cells in a contour to perform test (default 10).
#' @param p_adjust_method Method for multiple testing correction (default "fdr").
#'
#' @return A data.frame with columns: contour, gene/signature, fold_change, p_value, p_adj.
#' @export
#'
#' @examples
#' enrichment_df <- enrichment_in_contours(seurat_obj,
#'                                       genes = c("COL1A1", "ACTA2"),
#'                                       contour_col = "assigned_region")
enrichment_in_contours <- function(seurat_obj,
                                   genes,
                                   contour_col = "assigned_region",
                                   assay = "RNA",
                                   slot = "data",
                                   test_method = c("wilcox", "t.test"),
                                   min_cells = 10,
                                   p_adjust_method = "fdr") {
  test_method <- match.arg(test_method)
  meta <- seurat_obj@meta.data
  
  if (!contour_col %in% colnames(meta)) {
    stop(paste("Contour column", contour_col, "not found in metadata."))
  }
  
  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Prepare expression matrix or signature scores
  if (is.list(genes)) {
    # Calculate average expression per gene set (signature)
    expr_mat <- sapply(names(genes), function(sig_name) {
      sig_genes <- genes[[sig_name]]
      valid_genes <- intersect(sig_genes, rownames(expr_data))
      if (length(valid_genes) == 0) return(rep(NA_real_, ncol(seurat_obj)))
      Matrix::colMeans(expr_data[valid_genes, , drop = FALSE])
    })
  } else {
    valid_genes <- intersect(genes, rownames(expr_data))
    if (length(valid_genes) == 0) stop("None of the specified genes found in expression data.")
    expr_mat <- as.matrix(expr_data[valid_genes, , drop = FALSE])
    if (length(valid_genes) == 1) {
      expr_mat <- matrix(expr_mat, nrow = 1,
                         dimnames = list(valid_genes, colnames(seurat_obj)))
    }
    expr_mat <- t(expr_mat)  # cells x genes
  }
  
  expr_df <- as.data.frame(expr_mat)
  expr_df[[contour_col]] <- meta[[contour_col]]
  
  contours <- unique(expr_df[[contour_col]])
  genes_or_sigs <- setdiff(colnames(expr_df), contour_col)
  
  results <- list()
  
  for (contour in contours) {
    in_cells <- expr_df[[contour_col]] == contour
    out_cells <- !in_cells
    
    # Skip if not enough cells
    if (sum(in_cells) < min_cells) next
    
    for (g in genes_or_sigs) {
      in_expr <- expr_df[in_cells, g]
      out_expr <- expr_df[out_cells, g]
      
      # Skip if all NA
      if (all(is.na(in_expr)) || all(is.na(out_expr))) next
      
      # Compute fold change: mean(in) / mean(out), avoid div by zero
      mean_in <- mean(in_expr, na.rm = TRUE)
      mean_out <- mean(out_expr, na.rm = TRUE)
      fold_change <- ifelse(mean_out == 0, NA_real_, mean_in / mean_out)
      
      # Perform test
      p_val <- tryCatch({
        if (test_method == "wilcox") {
          suppressWarnings(stats::wilcox.test(in_expr, out_expr)$p.value)
        } else {
          suppressWarnings(stats::t.test(in_expr, out_expr)$p.value)
        }
      }, error = function(e) NA_real_)
      
      results[[length(results) + 1]] <- data.frame(
        contour = contour,
        gene_or_signature = g,
        fold_change = fold_change,
        p_value = p_val,
        stringsAsFactors = FALSE
      )
    }
  }
  
  results_df <- do.call(rbind, results)
  results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  
  return(results_df)
}
