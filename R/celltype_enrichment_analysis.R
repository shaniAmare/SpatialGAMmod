#' Perform cell type enrichment analysis across spatial regions
#'
#' This function tests whether certain cell types are significantly enriched in specific regions
#' compared to all other regions, using Fisher's exact test (default) or Chi-squared test.
#'
#' @param seurat_obj A Seurat object containing region and cell type metadata.
#' @param region_col Metadata column defining spatial regions (e.g., "assigned_region").
#' @param celltype_col Metadata column defining cell types.
#' @param method Statistical test to use: "fisher" (default) or "chisq".
#' @param min_cells Minimum number of cells per group to consider for testing. Default is 5.
#' @param p_adjust_method Method for multiple testing correction (passed to `p.adjust`). Default: "fdr".
#'
#' @return A data.frame with columns: region, celltype, p_value, p_adj, odds_ratio, enrichment (TRUE/FALSE).
#' @export
#'
#' @examples
#' enrichment_df <- celltype_enrichment_analysis(seurat_obj,
#'                        region_col = "assigned_region",
#'                        celltype_col = "CellType")
celltype_enrichment_analysis <- function(seurat_obj,
                                         region_col = "assigned_region",
                                         celltype_col = "CellType",
                                         method = c("fisher", "chisq"),
                                         min_cells = 5,
                                         p_adjust_method = "fdr") {
  method <- match.arg(method)
  meta <- seurat_obj@meta.data
  
  if (!region_col %in% colnames(meta)) {
    stop(paste("Region column", region_col, "not found in metadata."))
  }
  
  if (!celltype_col %in% colnames(meta)) {
    stop(paste("Cell type column", celltype_col, "not found in metadata."))
  }
  
  region_levels <- unique(meta[[region_col]])
  celltype_levels <- unique(meta[[celltype_col]])
  
  results <- list()
  
  for (region in region_levels) {
    in_region <- meta[[region_col]] == region
    
    for (celltype in celltype_levels) {
      is_celltype <- meta[[celltype_col]] == celltype
      
      # Build contingency table:
      #                 In Region   Not In Region
      # Is Celltype         A             B
      # Not Celltype        C             D
      
      A <- sum(in_region & is_celltype)
      B <- sum(!in_region & is_celltype)
      C <- sum(in_region & !is_celltype)
      D <- sum(!in_region & !is_celltype)
      
      total_cells <- A + B + C + D
      
      # Skip if too small
      if ((A + C) < min_cells || (B + D) < min_cells) next
      
      contingency <- matrix(c(A, B, C, D), nrow = 2,
                            dimnames = list(
                              Region = c("In", "Out"),
                              CellType = c("Yes", "No")
                            ))
      
      test_result <- if (method == "fisher") {
        suppressWarnings(fisher.test(contingency))
      } else {
        suppressWarnings(chisq.test(contingency))
      }
      
      results[[length(results) + 1]] <- data.frame(
        region = region,
        celltype = celltype,
        p_value = test_result$p.value,
        odds_ratio = ifelse(method == "fisher", test_result$estimate, NA_real_),
        enrichment = A / (A + C) > B / (B + D)
      )
    }
  }
  
  results_df <- do.call(rbind, results)
  
  results_df$p_adj <- p.adjust(results_df$p_value, method = p_adjust_method)
  
  return(results_df)
}