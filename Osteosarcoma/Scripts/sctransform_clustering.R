library(Seurat)
library(ggplot2)
library(dplyr)

#' Process single sample with SCTransform and clustering
#' @param object Seurat object
#' @param sample_id Sample identifier
#' @param ncells Number of cells for sketching
#' @return Processed Seurat object
process_sample <- function(object, sample_id, ncells = 50000) {
  tryCatch({
    # Set assay and prepare data
    DefaultAssay(object) <- "Spatial.008um"
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    
    # Sketch data
    object <- SketchData(
      object = object,
      ncells = ncells,
      method = "LeverageScore",
      sketched.assay = "sketch"
    )
    
    # Switch to sketch assay and perform clustering
    DefaultAssay(object) <- "sketch"
    object <- FindVariableFeatures(object)
    object <- ScaleData(object)
    object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
    object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
    object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 0.5)
    object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", 
                      return.model = TRUE, dims = 1:50)
    
    # Project data
    object <- ProjectData(
      object = object,
      assay = "Spatial.008um",
      full.reduction = "full.pca.sketch",
      sketched.assay = "sketch",
      sketched.reduction = "pca.sketch",
      umap.model = "umap.sketch",
      dims = 1:50,
      refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
    )
    
    return(object)
  }, error = function(e) {
    message(sprintf("Error processing sample %s: %s", sample_id, e$message))
    return(NULL)
  })
}

#' Generate standard visualizations
#' @param object Processed Seurat object
#' @param sample_id Sample identifier
#' @return List of plots
generate_visualizations <- function(object, sample_id) {
  # Sketched clustering plot
  DefaultAssay(object) <- "sketch"
  Idents(object) <- "seurat_cluster.sketched"
  p1 <- DimPlot(object, reduction = "umap.sketch", label = FALSE) + 
    ggtitle(paste0(sample_id, " - Sketched clustering")) + 
    theme(legend.position = "bottom")
  
  # Projected clustering plot
  DefaultAssay(object) <- "Spatial.008um"
  Idents(object) <- "seurat_cluster.projected"
  p2 <- DimPlot(object, reduction = "full.umap.sketch", label = FALSE) + 
    ggtitle(paste0(sample_id, " - Projected clustering")) + 
    theme(legend.position = "bottom")
  
  # Spatial dimension plot
  p3 <- SpatialDimPlot(object, label = TRUE, repel = TRUE, label.size = 4)
  
  return(list(sketch = p1, projected = p2, spatial = p3))
}

# Main processing pipeline
main <- function() {
  # Load all RDS files from directory
  files <- list.files("rds_files", pattern = "\\.rds$", full.names = TRUE)
  processed_objects <- list()
  
  for (file in files) {
    sample_id <- tools::file_path_sans_ext(basename(file))
    message(sprintf("Processing %s...", sample_id))
    
    # Load and process
    object <- readRDS(file)
    processed <- process_sample(object, sample_id)
    
    if (!is.null(processed)) {
      # Generate and save visualizations
      plots <- generate_visualizations(processed, sample_id)
      ggsave(
        sprintf("plots/%s_clustering.pdf", sample_id),
        plot = plots$sketch | plots$projected,
        width = 12, height = 6
      )
      ggsave(
        sprintf("plots/%s_spatial.pdf", sample_id),
        plot = plots$spatial,
        width = 8, height = 8
      )
      
      # Store processed object
      processed_objects[[sample_id]] <- processed
      saveRDS(processed, sprintf("processed/%s_processed.rds", sample_id))
    }
  }
  
  return(processed_objects)
}

# Create output directories if they don't exist
dir.create("plots", showWarnings = FALSE)
dir.create("processed", showWarnings = FALSE)

# Run pipeline
processed_objects <- main()