library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Configuration
SAMPLE_DIRS <- list(
  JY4 = "/Users/cameronwalker/Documents/BCM/VisiumHD/processed_data/VisiumHD_01_D1_JY4/outs",
  JY8 = "/Users/cameronwalker/Documents/BCM/VisiumHD/processed_data/JY8/outs",
  JY10 = "/Users/cameronwalker/Documents/BCM/visiumHD/processed_data/VisiumHD_01_A1_JY10/outs",
  JY23 = "/Users/cameronwalker/Documents/BCM/VisiumHD/processed_data/JY23/outs"
)

#' Load and initialize Visium data
#' @param data_dir Path to 10X output directory
#' @param bin_sizes Vector of bin sizes
#' @return Initialized Seurat object
load_visium_data <- function(data_dir, bin_sizes = c(8, 16)) {
  Load10X_Spatial(data_dir = data_dir, bin.size = bin_sizes)
}

#' Create QC plots for a Seurat object
#' @param object Seurat object
#' @param sample_name Name of the sample for plot titles
#' @return List containing violin and spatial plots
create_qc_plots <- function(object, sample_name) {
  DefaultAssay(object) <- "Spatial.008um"
  
  vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
    theme(axis.text = element_text(size = 4)) + 
    NoLegend() +
    ggtitle(paste0(sample_name, " - Count Distribution"))
  
  count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
    theme(legend.position = "right") +
    ggtitle(paste0(sample_name, " - Spatial Distribution"))
  
  return(list(violin = vln.plot, spatial = count.plot))
}

#' Normalize data for both 8um and 16um bins
#' @param object Seurat object
#' @return Normalized Seurat object
normalize_binned_data <- function(object) {
  assays <- c("Spatial.008um", "Spatial.016um")
  for (assay in assays) {
    DefaultAssay(object) <- assay
    object <- NormalizeData(object)
  }
  return(object)
}

# Main processing pipeline
process_samples <- function() {
  seurat_objects <- list()
  
  for (sample_name in names(SAMPLE_DIRS)) {
    # Load data
    message(sprintf("Processing sample %s...", sample_name))
    object <- load_visium_data(SAMPLE_DIRS[[sample_name]])
    
    # Generate QC plots
    qc_plots <- create_qc_plots(object, sample_name)
    print(qc_plots$violin | qc_plots$spatial)
    
    # Normalize data
    object <- normalize_binned_data(object)
    
    # Store processed object
    seurat_objects[[sample_name]] <- object
  }
  
  return(seurat_objects)
}

# Execute pipeline
seurat_objects <- process_samples()

# Save processed objects
saveRDS(seurat_objects, "processed_seurat_objects.rds")