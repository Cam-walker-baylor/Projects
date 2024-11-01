# Load necessary libraries
library(Seurat)      # For single-cell analysis
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(patchwork)   # For combining plots
library(future)      # For parallelization

# Set up parallel processing
plan("multicore", workers = 4)  # Adjust 'workers' based on available cores
options(future.globals.maxSize = 50 * 1024^3)  # Set max memory (50 GB in this example)

# Define the base directory where your data is stored
base_dir <- "/home/walker/data/projects/visiumHD/scRef/data"

# List of sample identifiers
sample_ids <- c("BC10", "BC11", "BC16", "BC17", "BC2", "BC20", "BC21", "BC22", "BC3", "BC5", "BC6")

# Initialize an empty list to store Seurat objects
seurat_list <- list()

# Loop over each sample to read data and create Seurat objects
for (sample_id in sample_ids) {
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

# Merge the Seurat objects into one integrated object
combined_seurat <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_ids)

# Perform QC (adjust thresholds as needed)
min_features <- 200
max_features <- 6000
max_mito <- 10  # Percent mitochondrial genes

# Subset the combined Seurat object based on QC metrics
combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mito)

# Normalize the data using SCTransform
combined_seurat <- SCTransform(combined_seurat, vars.to.regress = "percent.mt", verbose = FALSE)

# Perform PCA for dimensionality reduction
combined_seurat <- RunPCA(combined_seurat, npcs = 50)

elbow_plot <- ElbowPlot(combined_seurat, ndims = 50)
ggsave("elbow_plot.png", plot = elbow_plot, width = 8, height = 6, units = "in", dpi = 300)

# Find neighbors
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)

# Find clusters (adjust resolution as needed)

# Save the current plan
old_plan <- plan()

# Switch to sequential processing
plan("sequential")

# Run clustering
combined_seurat <- FindClusters(combined_seurat, resolution = 0.1)

# Restore the original plan
plan(old_plan)

# Run UMAP for further dimensionality reduction
combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)

# Create UMAP plot and store it in a variable
umap_plot <- DimPlot(combined_seurat, 
                     reduction = "umap",
                     label = TRUE,           # show cluster labels
                     label.size = 4,         # size of labels
                     pt.size = 0.5)          # size of points

# Save the plot
ggsave("umap_clusters.png", 
       plot = umap_plot, 
       width = 10, 
       height = 8, 
       dpi = 300, 
       units = "in")

# First, prepare the SCT data for marker analysis
combined_seurat <- PrepSCTFindMarkers(combined_seurat)

# Then run FindAllMarkers
all_markers <- FindAllMarkers(combined_seurat,
                              only.pos = TRUE,          
                              min.pct = 0.25,           
                              logfc.threshold = 1.00,    
                              assay = "SCT")            # specify SCT assay

# Sort and save results
all_markers <- all_markers %>%
  arrange(cluster, p_val_adj)

# Save to CSV
write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)













# Define gene sets (cell type markers)
Osteoblastic_OS_cells <- c("COL1A1", "CDH11", "RUNX2")
Proliferating_osteoblastic_OS_cells <- c("TOP2A", "PCNA", "MKI67")
Chondroblastic_OS_cells <- c("ACAN", "COL2A1", "SOX9")
Osteoclastic_cells <- c("CTSK", "MMP9")
TILs_T_and_NK_cells <- c("IL7R", "CD3D", "NKG7")
Myeloid_cells <- c("CD74", "CD14", "FCGR3A")
Fibroblasts <- c("COL1A1", "LUM", "DCN")
Pericytes <- c("ACTA2", "RGS5")
MSCs <- c("CXCL12", "SFRP2", "MME")
Myoblasts <- c("MYLPF", "MYL1")
Endothelial_cells <- c("PECAM1", "VWF")

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

# Add module scores for each gene set to the combined Seurat object
combined_seurat <- AddModuleScore(combined_seurat, features = gene_sets, name = names(gene_sets))

# Save the combined Seurat object with module scores
saveRDS(combined_seurat, file = "/data/u/walker/projects/osteosarcoma/data/human/combined_seurat_with_module_scores.rds")

# Extract module scores and cluster information
module_scores <- FetchData(combined_seurat, vars = paste0(names(gene_sets), "1"))
module_scores$cell_barcode <- rownames(module_scores)
module_scores$cluster <- Idents(combined_seurat)

# Save module scores to CSV
write.csv(module_scores, file = "/data/u/walker/projects/osteosarcoma/data/human/module_scores.csv", row.names = FALSE)

# Visualize module scores on UMAP
FeaturePlot(combined_seurat, features = paste0(names(gene_sets), "1"), reduction = "umap", ncol = 4)





######################################
# Define gene sets (cell type markers)
gene_sets <- list(
  Osteoblastic_OS_cells = c("COL1A1", "CDH11", "RUNX2"),
  Proliferating_osteoblastic_OS_cells = c("TOP2A", "PCNA", "MKI67"),
  Chondroblastic_OS_cells = c("ACAN", "COL2A1", "SOX9"),
  Osteoclastic_cells = c("CTSK", "MMP9"),
  TILs_T_and_NK_cells = c("IL7R", "CD3D", "NKG7"),
  Myeloid_cells = c("CD74", "CD14", "FCGR3A"),
  Fibroblasts = c("COL1A1", "LUM", "DCN"),
  Pericytes = c("ACTA2", "RGS5"),
  MSCs = c("CXCL12", "SFRP2", "MME"),
  Myoblasts = c("MYLPF", "MYL1"),
  Endothelial_cells = c("PECAM1", "VWF")
)

# Add module scores for each gene set to the combined Seurat object
combined_seurat <- AddModuleScore(combined_seurat, features = gene_sets, name = names(gene_sets))

# Verify that module scores were added
print(head(combined_seurat[[]]))

# Extract module scores and cluster information
module_score_names <- paste0(names(gene_sets), "1")
module_scores <- FetchData(combined_seurat, vars = c(module_score_names, "seurat_clusters"))
module_scores$cell_barcode <- rownames(module_scores)

# Save module scores to CSV
write.csv(module_scores, file = "module_scores.csv", row.names = FALSE)

# Visualize module scores on UMAP
feature_plots <- FeaturePlot(combined_seurat, features = module_score_names, reduction = "umap", ncol = 3)

# Print the feature plots
print(feature_plots)

# Save the feature plots
ggsave("module_scores_umap.pdf", feature_plots, width = 20, height = 30, units = "in")

# Save the combined Seurat object with module scores
saveRDS(combined_seurat, file = "combined_seurat_with_module_scores.rds")

# Print a summary of the Seurat object
print(combined_seurat)