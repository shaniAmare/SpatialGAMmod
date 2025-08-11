# SpatialGAMmod


**SpatialGAMmod** is an R package designed for spatial transcriptomics data analysis. It provides tools to model spatial gene expression gradients using Generalized Additive Models (GAMs), perform spatial autocorrelation analyses, compare tumor vs non-tumor regions, analyze cell type proportions, and visualize spatial patterns effectively.

*NOTE: This is only a side project of mine, so I might be slow to respond to any issues.*


---

## Features

- Fit spatial GAMs to model gene or signature expression gradients.
- Calculate spatial autocorrelation statistics (Moran's I).
- Analyze spatial correlation decay with distance.
- Compute and visualize cell type proportions and enrichments.
- Compare tumor vs non-tumor contour regions for differential expression.
- Run pseudotime inference within spatial regions.
- Export cells by spatial region and perform various contour-based analyses.
- Comprehensive visualization functions including heatmaps, barplots, density maps, and contour plots.

---

## Installation

### 1. Install dependencies

The package relies on several R packages. You can install the required dependencies from CRAN and Bioconductor as follows:

```r
# CRAN packages
install.packages(c(
  "Seurat",
  "mgcv",
  "ggplot2",
  "viridis",
  "dplyr",
  "tidyr",
  "pheatmap",
  "Matrix",
  "SingleCellExperiment",
  "slingshot",
  "spdep"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "slingshot", "Matrix"))

# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("shaniAmare/SpatialGAMmod")
```

### 2. Basic Usage

#### Fit spatial GAM gradients for genes

```r
library(SpatialGAMmod)

result <- spatial_gam_gradient(
  seurat_obj = your_seurat_obj,
  genes = c("ACTA2", "COL1A1"),
  split_by = "tumor_core_location",
  plot_results = TRUE
)

# View GAM gradient plots
print(result$plots[[1]])
```

#### Calculate spatial autocorrelation (Moran’s I)

```r
moran_i <- spatial_autocorrelation(
  seurat_obj = your_seurat_obj,
  gene = "ACTA2",
  coord_cols = c("imagecol", "imagerow"),
  neighbors_k = 6
)
print(moran_i)
```

#### Compare tumor vs non-tumor regions

```r
diff_res <- tumour_vs_non_tumour_contour(
  seurat_obj = your_seurat_obj,
  contour_col = "assigned_region",
  tumour_regions = c("tumour_core"),
  non_tumour_regions = c("stroma", "immune"),
  genes = c("ACTA2", "COL1A1")
)

# View volcano plot
print(diff_res$plot)
```

#### Run pseudotime within regions

```r
seurat_obj <- run_pseudotime_within_region(
  seurat_obj = your_seurat_obj,
  region_col = "assigned_region",
  reduction = "pca",
  dims = 1:10,
  cluster_col = "celltype"
)
```

### Complete User Manual

For a detailed user manual and workflow guide, see the [User Manual](vignettes/user_manual.html).

### Package Functions Overview
```css
- aggregate_expression_per_region()
- assign_contour_regions()
- calculate_celltype_proportions()
- celltype_distribution_plots()
- celltype_enrichment_analysis()
- compare_expression_across_contours()
- compute_diversity_index()
- compute_entropy_per_region()
- contour_variability_by_group()
- enrichment_in_contours()
- export_cells_by_region()
- fit_spatial_gam()
- gene_signature_by_contour()
- localised_interactions_summary()
- map_contour_to_cells()
- marker_density_maps()
- plot_celltype_barplot()
- plot_celltype_heatmap()
- plot_expression_density_per_region()
- plot_gam_contours()
- predict_gam_to_grid()
- proportion_heatmaps()
- run_pseudotime_within_region()
- spatial_autocorrelation()
- spatial_correlation_with_distance()
- spatial_gam_gradient()
- tumour_vs_non_tumour_contour()
```
---

### Citation
Please cite this package if used in published research.

### Contact
For issues, feature requests, or questions, please open an issue on GitHub or contact shani.amarasinghe@monash.edu.

*Happy spatial analysis!*
