#' Calculate and summarize gene signature scores by spatial contour region
#'
#' Computes average expression-based gene signature scores per cell and summarizes
#' these scores across spatial contour regions.
#'
#' @param seurat_obj A Seurat object with expression data and metadata.
#' @param gene_signatures Named list of gene vectors, each representing a signature.
#' @param contour_col Metadata column defining spatial contours (e.g., "assigned_region").
#' @param assay Assay to use for expression data (default "RNA").
#' @param slot Slot in assay to extract expression from (default "data").
#' @param summary_fun Function to summarize signature scores per contour (default mean).
#' @param min_genes Minimum number of signature genes required per signature per cell (default 3).
#' @param return_per_cell Logical; if TRUE, return per-cell signature scores as well (default FALSE).
#'
#' @return A list with elements:
#'   - contour_summary: data.frame summarizing signature scores per contour (rows) and signature (columns).
#'   - per_cell_scores: data.frame of per-cell signature scores with metadata (if return_per_cell=TRUE).
#' @export
#'
#' @examples
#' signatures <- list(
#'   Fibroblast = c("COL1A1", "ACTA2", "FAP"),
#'   Immune = c("CD3D", "CD79A", "PTPRC")
#' )
#' result <- gene_signature_by_contour(seurat_obj,
#'                                    gene_signatures = signatures,
#'                                    contour_col = "assigned_region",
#'                                    return_per_cell = TRUE)
gene_signature_by_contour <- function(seurat_obj,
                                      gene_signatures,
                                      contour_col = "assigned_region",
                                      assay = "RNA",
                                      slot = "data",
                                      summary_fun = mean,
                                      min_genes = 3,
                                      return_per_cell = FALSE) {
  meta <- seurat_obj@meta.data
  
  # Check contour_col presence
  if (!contour_col %in% colnames(meta)) {
    stop(paste0("Contour column '", contour_col, "' not found in Seurat metadata."))
  }
  
  # Validate gene_signatures list
  if (!is.list(gene_signatures) || is.null(names(gene_signatures))) {
    stop("gene_signatures must be a named list of character vectors (genes).")
  }
  
  # Extract expression matrix (genes x cells)
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  # Prepare container for per-cell scores
  per_cell_scores <- matrix(NA_real_, nrow = ncol(seurat_obj), ncol = length(gene_signatures))
  colnames(per_cell_scores) <- names(gene_signatures)
  rownames(per_cell_scores) <- colnames(seurat_obj)
  
  # Calculate per-cell signature scores
  for (sig_name in names(gene_signatures)) {
    sig_genes <- gene_signatures[[sig_name]]
    valid_genes <- intersect(sig_genes, rownames(expr_mat))
    
    if (length(valid_genes) < min_genes) {
      warning(sprintf("Signature '%s' has fewer than %d genes found in data. Skipping.", sig_name, min_genes))
      next
    }
    
    # Calculate mean expression per cell for signature genes
    sig_expr <- Matrix::colMeans(expr_mat[valid_genes, , drop = FALSE])
    per_cell_scores[, sig_name] <- sig_expr
  }
  
  # Build data.frame for per-cell scores and metadata
  per_cell_scores_df <- data.frame(per_cell_scores, check.names = FALSE)
  per_cell_scores_df[[contour_col]] <- meta[[contour_col]]
  per_cell_scores_df$cell <- rownames(per_cell_scores_df)
  
  # Summarize per contour region
  contour_summary <- per_cell_scores_df |>
    tidyr::pivot_longer(
      cols = names(gene_signatures),
      names_to = "signature",
      values_to = "score"
    ) |>
    dplyr::group_by(!!rlang::sym(contour_col), signature) |>
    dplyr::summarise(
      summary_score = summary_fun(score, na.rm = TRUE),
      n_cells = sum(!is.na(score)),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = signature,
      values_from = summary_score
    ) |>
    dplyr::rename(region = !!rlang::sym(contour_col))
  
  if (return_per_cell) {
    return(list(
      contour_summary = contour_summary,
      per_cell_scores = per_cell_scores_df
    ))
  } else {
    return(list(
      contour_summary = contour_summary
    ))
  }
}