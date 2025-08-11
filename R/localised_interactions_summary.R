#' Summarize localized cell–cell interactions within spatial contours
#'
#' Given a data.frame of cell–cell interaction pairs and contour assignments,
#' this function calculates interaction counts or frequencies between cell types,
#' summarized per contour region.
#'
#' @param interactions_df A data.frame with columns:
#'   - cell1: first cell barcode/ID
#'   - cell2: interacting neighbor cell barcode/ID
#'   - celltype1: cell type or cluster label for cell1
#'   - celltype2: cell type or cluster label for cell2
#'   - contour: contour/region assignment of cell1 (or interaction location)
#' @param region_col Name of the column in `interactions_df` containing contour/region assignment (default "contour").
#' @param interaction_type Optional character; type of interaction to filter by (if available).
#' @param summarize_by One of "counts" (default) or "frequency" (proportion of interactions).
#'
#' @return A data.frame with rows as contour regions and columns as interaction pairs (e.g., "Fibroblast->Immune"),
#'         containing counts or frequencies of interactions per region.
#' @export
#'
#' @examples
#' # Sample input format: interactions_df with columns cell1, cell2, celltype1, celltype2, contour
#' interaction_summary <- localised_interactions_summary(interactions_df)
localised_interactions_summary <- function(interactions_df,
                                           region_col = "contour",
                                           interaction_type = NULL,
                                           summarize_by = c("counts", "frequency")) {
  summarize_by <- match.arg(summarize_by)
  
  required_cols <- c("cell1", "cell2", "celltype1", "celltype2", region_col)
  if (!all(required_cols %in% colnames(interactions_df))) {
    stop("interactions_df must contain columns: cell1, cell2, celltype1, celltype2, and the region_col.")
  }
  
  df <- interactions_df
  
  # Optionally filter by interaction_type column if provided
  if (!is.null(interaction_type)) {
    if (!"interaction_type" %in% colnames(df)) {
      warning("interaction_type specified but not found in interactions_df; ignoring filter.")
    } else {
      df <- df[df$interaction_type == interaction_type, , drop = FALSE]
    }
  }
  
  # Create interaction pair labels (e.g., "Fibroblast->Immune")
  df$interaction_pair <- paste(df$celltype1, "->", df$celltype2)
  
  # Count interactions per region and interaction pair
  counts_df <- df |>
    dplyr::group_by(!!rlang::sym(region_col), interaction_pair) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  # Spread to wide format: rows = region, columns = interaction pairs
  counts_wide <- tidyr::pivot_wider(counts_df,
                                    names_from = interaction_pair,
                                    values_from = count,
                                    values_fill = 0)
  
  if (summarize_by == "frequency") {
    # Calculate total interactions per region for normalization
    total_interactions <- rowSums(counts_wide[ , -1, drop = FALSE])
    counts_wide[ , -1] <- sweep(counts_wide[ , -1, drop = FALSE], 1, total_interactions, FUN = "/")
  }
  
  # Rename first column to "region" for clarity
  colnames(counts_wide)[1] <- "region"
  
  return(counts_wide)
}