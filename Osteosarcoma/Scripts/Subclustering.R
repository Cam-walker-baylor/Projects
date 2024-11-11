library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)
library(tidyr)
library(igraph)
library(reshape2)
library(SpatialExperiment)
library(scater)
library(SeuratDisk)
library(stringr)  # For string manipulation
library(SingleCellExperiment)

# Clear the environment to avoid variable conflicts
rm(list = ls())

# **Step 1: Load Data and Prepare Seurat Object**

# Load the Seurat object for JY10
JY10 <- readRDS("rds_files/JY10_object.rds")

# Print the Seurat object summary
print(JY10)

# Get the barcodes (cell names) from the Seurat object
JY10_barcodes <- Cells(JY10)

# Display the first few barcodes
head(JY10_barcodes)

# Read in the cluster assignments CSV file for JY10
JY10_cluster_assignments <- read.csv("JY10_ANNOTATED.csv", header = TRUE, stringsAsFactors = FALSE)

# Display the first few rows of cluster assignments
head(JY10_cluster_assignments)

# First, let's make sure the barcodes in cluster assignments match those in the Seurat object
# and create a clean cell type annotation vector

# Create a named vector of cell type assignments
cell_type_assignments <- setNames(
  JY10_cluster_assignments$Graph.based,
  JY10_cluster_assignments$Barcode
)

# Get unique cell types (excluding empty/NA values)
unique_cell_types <- unique(cell_type_assignments[cell_type_assignments != ""])
unique_cell_types <- unique_cell_types[!is.na(unique_cell_types)]

# Create a list to store all Seurat objects
seurat_objects_by_celltype <- list()

# Function to create subset Seurat object for each cell type
create_celltype_objects <- function(seurat_obj, cell_assignments, cell_type) {
  # Get barcodes for the specific cell type
  cell_barcodes <- names(cell_assignments[cell_assignments == cell_type])
  
  # Create subset
  subset_obj <- subset(seurat_obj, cells = cell_barcodes)
  
  # Add metadata about the cell type
  subset_obj[["cell_type"]] <- cell_type
  
  return(subset_obj)
}

# Create separate Seurat objects for each cell type
for (cell_type in unique_cell_types) {
  # Create clean name for object (replace spaces and special characters)
  clean_name <- gsub("[^[:alnum:]]", "_", cell_type)
  
  # Create the subset object
  seurat_objects_by_celltype[[clean_name]] <- create_celltype_objects(
    JY10, 
    cell_type_assignments, 
    cell_type
  )
  
  # Print summary for this cell type
  cat(sprintf("\nCreated Seurat object for %s:", cell_type))
  print(seurat_objects_by_celltype[[clean_name]])
}

# Optional: Save individual Seurat objects
save_objects <- function(seurat_objects_list, output_dir = "seurat_objects_by_celltype") {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save each object
  for (obj_name in names(seurat_objects_list)) {
    output_file <- file.path(output_dir, paste0(obj_name, ".rds"))
    saveRDS(seurat_objects_list[[obj_name]], output_file)
    cat(sprintf("\nSaved %s to %s", obj_name, output_file))
  }
}

# Uncomment the following line to save the objects
save_objects(seurat_objects_by_celltype)

# Function to process a spatial Seurat subset directly
process_spatial_subset <- function(seurat_obj, 
                                   dims.use = 1:30,
                                   resolution = 0.1) {
  
  # Set default assay to Spatial.008um
  DefaultAssay(seurat_obj) <- "Spatial.008um"
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  # Scale only the variable features
  var_features <- VariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = var_features)
  
  # Run PCA on variable features
  seurat_obj <- RunPCA(seurat_obj, features = var_features)
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = dims.use)
  
  # Find neighbors
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims.use)
  
  # Find clusters
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  return(seurat_obj)
}

# Function to create visualizations
create_spatial_plots <- function(seurat_obj, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # UMAP clustering plot
  p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle("UMAP clustering") +
    theme(legend.position = "right")
  
  # Save UMAP plot
  pdf(file.path(output_dir, "umap_clusters.pdf"), width = 10, height = 8)
  print(p1)
  dev.off()
  
  # Spatial dimension plot
  pdf(file.path(output_dir, "spatial_clusters.pdf"), width = 10, height = 8)
  print(SpatialDimPlot(seurat_obj, label = TRUE, repel = TRUE, label.size = 4))
  dev.off()
  
  # Find markers
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  
  if(nrow(markers) > 0) {
    # Save markers
    write.csv(markers, file.path(output_dir, "cluster_markers.csv"), row.names = FALSE)
    
    # Create and save heatmap of top markers
    markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = 5) %>%
      ungroup() -> top5
    
    # Scale only the genes we need for the heatmap
    seurat_obj <- ScaleData(seurat_obj, 
                            features = unique(top5$gene), 
                            assay = "Spatial.008um")
    
    p_markers <- DoHeatmap(seurat_obj, 
                           features = unique(top5$gene), 
                           size = 2.5) +
      theme(axis.text = element_text(size = 5.5)) +
      NoLegend()
    
    pdf(file.path(output_dir, "markers_heatmap.pdf"), width = 12, height = 8)
    print(p_markers)
    dev.off()
    
    # Create spatial feature plots for top markers
    pdf(file.path(output_dir, "spatial_features.pdf"), width = 12, height = 8)
    for(gene in unique(top5$gene)) {
      print(SpatialFeaturePlot(seurat_obj, features = gene) + 
              theme(legend.position = "right"))
    }
    dev.off()
  }
}

########### Process single Object
# Test with one object
test_name <- "Myeloid_cells"
test_obj <- seurat_objects_by_celltype[[test_name]]
processed_test <- process_spatial_subset(test_obj)
create_spatial_plots(processed_test, "test_output")


# Process all objects
for (obj_name in names(seurat_objects_by_celltype)) {
  cat(sprintf("\nProcessing %s...\n", obj_name))
  output_dir <- file.path("analysis_results", obj_name)
  
  tryCatch({
    # Get current object
    current_obj <- seurat_objects_by_celltype[[obj_name]]
    
    # Process object
    processed_obj <- process_spatial_subset(current_obj)
    
    # Create plots
    create_spatial_plots(processed_obj, output_dir)
    
    # Save processed object
    saveRDS(processed_obj, file.path(output_dir, paste0(obj_name, "_processed.rds")))
    
    # Print basic stats
    cat(sprintf("Completed processing %s:\n", obj_name))
    cat(sprintf("Number of cells: %d\n", ncol(processed_obj)))
    cat(sprintf("Number of clusters: %d\n", length(unique(Idents(processed_obj)))))
    
  }, error = function(e) {
    cat(sprintf("Error processing %s: %s\n", obj_name, e$message))
  })
}