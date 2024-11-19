# Load required libraries
library(Seurat)      # For single-cell analysis
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(patchwork)   # For combining plots
library(future)      # For parallelization

# Set up parallel processing
plan("multicore", workers = 4)  # Adjust 'workers' based on available cores
options(future.globals.maxSize = 50 * 1024^3)  # Set max memory (50 GB in this example)

# Define the base directory where your data is stored
base_dir <- "/data/u/walker/projects/visiumHD/scRef/data"

# List of sample identifiers
sample_ids <- c("BC10", "BC11", "BC16", "BC17", "BC2", "BC20", "BC21", "BC22", "BC3", "BC5", "BC6")

# Initialize an empty list to store Seurat objects
seurat_list <- list()

# Loop over each sample to read data and create Seurat objects
for (sample_id in sample_ids) {
  # Construct the path to the sample's data directory
  data_dir <- file.path(base_dir, sample_id)
  
  # Read the 10X Genomics data (Cell Ranger output)
  seurat_obj <- Read10X(data.dir = data_dir)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = sample_id, min.cells = 3, min.features = 200)
  
  # Add mitochondrial percentage to metadata
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Store the Seurat object in the list
  seurat_list[[sample_id]] <- seurat_obj
}

# Define QC thresholds
min_features <- 200
max_features <- 6000
max_mito <- 10  # Percent mitochondrial genes

for (sample_id in names(seurat_list)) {
  seurat_obj <- seurat_list[[sample_id]]
  
  # Visualize QC metrics before filtering
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Filter cells based on QC metrics
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mito)
  
  # Update the object in the list
  seurat_list[[sample_id]] <- seurat_obj
}

for (sample_id in names(seurat_list)) {
  seurat_obj <- seurat_list[[sample_id]]
  
  # Normalize with SCTransform
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  
  # Update the object in the list
  seurat_list[[sample_id]] <- seurat_obj
}

# Select features for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# Prepare objects for integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)

# Integrate data
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Set the default assay to integrated
DefaultAssay(seurat_integrated) <- "integrated"

# Run PCA
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)

# Examine PCA results
ElbowPlot(seurat_integrated)

# Run UMAP
seurat_integrated03 <- RunUMAP(seurat_integrated, dims = 1:8)

# Find neighbors
seurat_integrated03 <- FindNeighbors(seurat_integrated03, dims = 1:8)

# Find clusters
seurat_integrated03 <- FindClusters(seurat_integrated03, resolution = 1.0)

# UMAP plot colored by cluster
DimPlot(seurat_integrated03, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

ggsave("dimplot_umap.png", 
       plot = DimPlot(seurat_integrated03, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend(),
       width = 10, height = 8, units = "in", dpi = 300)

saveRDS(seurat_integrated03, file = "seurat_integrated_clusters03.rds")

# Find markers for all clusters
markers04 <- FindAllMarkers(
  object = seurat_integrated03,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1.00
)

# Display top markers for each cluster
markers04 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Save markers to CSV
write.csv(markers04, file = "cluster_markers_dim8_res10.csv", row.names = FALSE)

# Example usage of AddModuleScore
# Define gene sets (replace with your own gene sets)
Osteoblastic_OS_cells <- c("Col1a1", "Cdh11", "Runx2")
Proliferating_osteoblastic_OS_cells <- c("Top2a", "Pcna", "Mki67")
Chondroblastic_OS_cells <- c("Acan", "Col2a1", "Sox9")
Osteoclastic_cells <- c("Ctsk", "Mmp9")
TILs_T_and_NK_cells <- c("Il7r", "Cd3d", "Nkg7")
Myeloid_cells <- c("Cd74", "Cd14", "Fcgr3a")
Fibroblasts <- c("Col1a1", "Lum", "Dcn")
Pericytes <- c("Acta2", "Rgs5")
MSCs <- c("Cxcl12", "Sfrp2", "Mme")
Myoblasts <- c("Mylpf", "Myl1")
Endothelial_cells <- c("Pecam1", "Vwf")

# Create a list of gene sets
gene_sets <- list(
  Osteoblastic_OS_cells = Osteoblastic_OS_cells,
  Proliferating_osteoblastic_OS_cells = Proliferating_osteoblastic_OS_cells,
  Chondroblastic_OS_cells = Chondroblastic_OS_cells,
  Osteoclastic_cells = Osteoclastic_cells,
  TILs_T_and_NK_cells = TILs_T_and_NK_cells,
  Myeloid_cells = Myeloid_cells,
  Fibroblasts = Fibroblasts,
  Pericytes = Pericytes,
  MSCs = MSCs,
  Myoblasts = Myoblasts,
  Endothelial_cells = Endothelial_cells
)

# Add module scores for each gene set
seurat_integrated04 <- AddModuleScore(seurat_integrated03, 
                                      features = gene_sets,
                                      name = names(gene_sets))
# Extract module scores from the Seurat object
module_scores <- FetchData(seurat_integrated03, 
                           vars = paste0(names(gene_sets), "1"))

# Add cell barcodes and cluster information
module_scores$cell_barcode <- rownames(module_scores)
module_scores$cluster <- Idents(seurat_integrated03)

# Save module scores to CSV
write.csv(module_scores, file = "module_scores.csv", row.names = FALSE)

# Visualize module scores on UMAP
FeaturePlot(seurat_integrated04, 
            features = paste0(names(gene_sets), "1"),
            reduction = "umap",
            ncol = 4)  # Adjust ncol as needed for better layout

# Save the updated object
saveRDS(seurat_integrated03, file = "seurat_integrated_clusters04_with__module_scores.rds")