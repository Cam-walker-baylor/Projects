# Visium Analysis Pipeline

A modular pipeline for analyzing Visium spatial transcriptomics data.

## Scripts Overview

### 1. Preprocessing Pipeline (`preprocess_visium.R`)
- Loads and initializes Visium data from Cell Ranger output
- Performs QC and normalization
- Generates QC visualizations
- Creates normalized Seurat objects ready for analysis

Usage:
```r
seurat_objects <- process_samples()
```

### 2. SCTransform and Clustering (`sctransform_clustering.R`)
- Processes normalized Seurat objects
- Performs data sketching for large datasets
- Executes clustering analysis
- Creates UMAP and spatial visualizations

Usage:
```r
processed_objects <- main()
```

### 3. Tumor Boundary Detection (`tumor_boundary.R`)
- Identifies tumor boundaries using KDE
- Processes hypoxia-labeled cells
- Creates contours and polygons
- Generates boundary visualizations

Usage:
```r
results <- run_tumor_boundary_pipeline(
  hypoxia_file = "data/hypoxia.csv",
  seurat_file = "data/seurat_obj.rds",
  output_prefix = "sample1"
)
```

### 4. Spatial Binning (`spatial_binning.R`)
- Performs spatial binning of Visium data
- Aggregates gene expression in bins
- Generates QC metrics and visualizations
- Creates binned Seurat objects

Usage:
```r
results <- run_binning_pipeline(samples, "visium_analysis")
```

### 5. Label Transfer (`label_transfer.R`)
- Transfers cell type labels between samples
- Handles OS and non-OS annotations separately
- Combines predictions with confidence scores
- Creates annotated Seurat objects

Usage:
```r
results <- run_label_transfer(
  reference_path = "ref.rds",
  query_path = "query.rds",
  annotations_path = "annotations.csv",
  output_dir = "output"
)
```

### 6. Colocalization Analysis (`colocalization.R`)
- Analyzes spatial relationships between cell types
- Performs statistical testing for colocalization
- Generates correlation matrices and heatmaps
- Creates network visualizations of cell type interactions

Usage:
```r
results <- run_colocalization_pipeline(
  sample_id = "sample1",
  panel_size = 500,
  overlap = 0.1,
  p_value_threshold = 0.00001
)
```

## Pipeline Dependencies
- Seurat (>= 4.0.0)
- tidyverse
- ggplot2
- MASS
- Matrix
- igraph
- dbscan

## Data Requirements
- Cell Ranger output directories
- Annotation files (CSV format)
- Hypoxia labeling data (for tumor boundary detection)

## Output Structure
```
output/
├── processed/
│   └── seurat_objects/
├── analysis/
│   ├── clustering/
│   ├── boundaries/
│   ├── binning/
│   └── colocalization/
├── plots/
└── metadata/
```

## Usage Notes
- Scripts can be run independently or as a pipeline
- Parameters configurable via function arguments
- Output directories created automatically
- Progress logging and error handling included