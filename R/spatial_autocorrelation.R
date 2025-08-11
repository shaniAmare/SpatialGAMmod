#' Calculate spatial autocorrelation (Moran's I) for gene expression in spatial data
#'
#' Computes Moran's I statistic for spatial autocorrelation of a numeric feature (e.g., gene expression)
#' using spatial neighbor weights derived from cell coordinates.
#'
#' @param seurat_obj A Seurat object with spatial coordinates in metadata.
#' @param gene Character; gene name to calculate autocorrelation for.
#' @param coord_cols Character vector of length 2 specifying spatial coordinate columns (default c("imagecol", "imagerow")).
#' @param assay Assay name to get expression from (default "RNA").
#' @param slot Slot in assay to use (default "data").
#' @param neighbors_k Integer; number of nearest neighbors for spatial weights (default 6).
#'
#' @return A numeric value of Moran's I statistic.
#' @importFrom spdep knearneigh knn2nb nb2listw moran.test
#' @export
#'
#' @examples
#' morans_i <- spatial_autocorrelation(seurat_obj, gene = "ACTA2")
spatial_autocorrelation <- function(seurat_obj,
                                    gene,
                                    coord_cols = c("imagecol", "imagerow"),
                                    assay = "RNA",
                                    slot = "data",
                                    neighbors_k = 6) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' needed for spatial autocorrelation calculations.")
  }
  
  meta <- seurat_obj@meta.data
  
  if (!all(coord_cols %in% colnames(meta))) {
    stop("Coordinate columns not found in Seurat metadata.")
  }
  
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste0("Gene '", gene, "' not found in Seurat object."))
  }
  
  # Extract coordinates and expression
  coords <- meta[, coord_cols]
  expr <- as.numeric(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, ])
  
  # Remove any NA cells
  valid <- complete.cases(coords) & !is.na(expr)
  coords <- coords[valid, , drop = FALSE]
  expr <- expr[valid]
  
  # Compute k-nearest neighbor graph
  library(spdep)
  
  knn <- spdep::knearneigh(coords, k = neighbors_k)
  nb <- spdep::knn2nb(knn)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Calculate Moran's I
  mi_test <- spdep::moran.test(expr, lw, zero.policy = TRUE)
  
  return(mi_test$estimate["Moran I statistic"])
}