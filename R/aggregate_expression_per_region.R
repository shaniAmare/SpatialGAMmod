#' Aggregate gene expression per spatial contour region
#'
#' This function computes the average expression of specified genes or gene signatures
#' within each spatial contour region assigned to cells in a Seurat object.
#'
#' @param seurat_obj A Seurat object containing spatial transcriptomics data with region labels.
#' @param genes Character vector of gene names to aggregate. Can also be a named list of gene sets.
#' @param region_col Name of the metadata column containing region or contour assignments.
#' @param assay Assay to pull expression data from. Defaults to "RNA".
#' @param slot Slot to extract expression from. Defaults to "data".
#' @param summary_fun Function to summarize expression per region (e.g., mean, median). Defaults to mean.
#' @param normalize Logical; if TRUE, expression is z-scored across regions per gene/gene set.
#'
#' @return A data.frame with aggregated expression values per region.
#' @export
#'
#' @examples
#' agg_df <- aggregate_expression_per_region(seurat_obj, genes = c("COL1A1", "ACTA2"), region_col = "contour_id")
aggregate_expression_per_region <- function(seurat_obj,
                                            genes,
                                            region_col = "contour_id",
                                            assay = "RNA",
                                            slot = "data",
                                            summary_fun = mean,
                                            normalize = FALSE) {
  # Check input
  if (!region_col %in% colnames(seurat_obj[[]])) {
    stop(paste("Column", region_col, "not found in metadata."))
  }
  
  expr_data <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # If genes is a named list (e.g., gene signatures), compute module scores
  if (is.list(genes)) {
    gene_expr <- sapply(names(genes), function(set_name) {
      genes_in_set <- genes[[set_name]]
      valid_genes <- intersect(genes_in_set, rownames(expr_data))
      if (length(valid_genes) == 0) return(rep(NA, ncol(expr_data)))
      colMeans(expr_data[valid_genes, , drop = FALSE])
    })
  } else {
    # Simple gene expression subset
    valid_genes <- intersect(genes, rownames(expr_data))
    if (length(valid_genes) == 0) stop("None of the specified genes are found in the expression matrix.")
    gene_expr <- expr_data[valid_genes, , drop = FALSE]
    if (length(valid_genes) == 1) {
      gene_expr <- matrix(gene_expr, nrow = 1,
                          dimnames = list(valid_genes, colnames(expr_data)))
    }
    gene_expr <- t(gene_expr)
  }
  
  gene_expr_df <- as.data.frame(gene_expr)
  gene_expr_df$region <- seurat_obj[[region_col]][, 1]
  
  # Aggregate
  agg_expr <- gene_expr_df |>
    dplyr::group_by(region) |>
    dplyr::summarise(across(where(is.numeric), summary_fun, na.rm = TRUE)) |>
    dplyr::ungroup()
  
  # Optionally normalize
  if (normalize) {
    numeric_cols <- names(agg_expr)[sapply(agg_expr, is.numeric)]
    agg_expr[numeric_cols] <- scale(agg_expr[numeric_cols])
  }
  
  return(as.data.frame(agg_expr))
}
