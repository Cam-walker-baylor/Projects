# Load Necessary Libraries
library(Seurat)         # Seurat V5
library(tidyverse)      # Includes readr, dplyr, ggplot2
library(dbscan)
library(Matrix)         # For sparse matrix operations
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

# Set options for debugging
options(stringsAsFactors = FALSE)

# Define a function to standardize barcodes
standardize_barcodes <- function(barcodes) {
  barcodes <- trimws(barcodes)
  barcodes <- toupper(barcodes)
  return(barcodes)
}

# Step 1: Load Data
cat("Loading tumor microenvironment labels...\n")
data <- read_csv("JY10_Tumor_Micro_ENV.csv")

# Adjust the column names if necessary
if (!all(c("Barcode", "2") %in% colnames(data))) {
  stop("The columns 'Barcode' and '2' must be present in JY10_Tumor_Micro_ENV.csv")
}

# Rename the '2' column to 'Label' for clarity
data <- data %>%
  rename(Label = `2`)

# Filter for barcodes labeled as "Tumor"
tumor_barcodes <- data %>%
  filter(Label == "Tumor") %>%
  select(Barcode)

# Standardize barcode formats
tumor_barcodes$Barcode <- standardize_barcodes(tumor_barcodes$Barcode)

# Step 2: Load Seurat Object and Extract Spatial Coordinates
cat("Loading Seurat object...\n")
JY10 <- readRDS("rds_files/JY10_object.rds")

# List available images
cat("Available images in the Seurat object:\n")
print(Images(JY10))

# Extract spatial coordinates for the desired image
image_name <- "slice1.008um"  # Adjust as necessary
if (!image_name %in% Images(JY10)) {
  stop(paste("Image", image_name, "not found in the Seurat object."))
}
spatial_coords <- GetTissueCoordinates(JY10, image = image_name)

# Standardize barcodes in spatial_coords
spatial_coords$cell <- standardize_barcodes(spatial_coords$cell)

# Step 3: Merge Tumor Labels with Spatial Coordinates
cat("Merging tumor labels with spatial coordinates...\n")
tumor_spatial_data <- spatial_coords %>%
  filter(cell %in% tumor_barcodes$Barcode) %>%
  mutate(Label = "Tumor")

# Debugging: Check number of matched tumor cells
cat("Number of tumor cells with spatial coordinates:", nrow(tumor_spatial_data), "\n")

# Step 3a: Generate a CSV with all barcodes labeled as "Tumor" or "Non-Tumor"
cat("Generating CSV with all barcodes labeled as 'Tumor' or 'Non-Tumor'...\n")
# Create a data frame with all barcodes from JY10
all_barcodes <- data.frame(Barcode = standardize_barcodes(Cells(JY10)))

# Label barcodes as "Tumor" or "Non-Tumor"
all_barcodes <- all_barcodes %>%
  mutate(Label = ifelse(Barcode %in% tumor_spatial_data$cell, "Tumor", "Non-Tumor"))

# Save to CSV
write_csv(all_barcodes, "JY10_all_barcodes_with_labels.csv")
cat("CSV file 'JY10_all_barcodes_with_labels.csv' generated.\n")

# Step 4: Visualize Initial Tumor-Labeled Cells
cat("Creating initial tumor plot...\n")
initial_tumor_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y)) +
  geom_point(color = "red", size = 0.5) +
  labs(title = "Spatial Distribution of Tumor-Labeled Cells") +
  theme_minimal()

ggsave("initial_tumor_plot.png", plot = initial_tumor_plot, width = 8, height = 6)
cat("Initial tumor plot saved as 'initial_tumor_plot.png'.\n")

# Step 5: Apply DBSCAN for Boundary Refinement
cat("Applying DBSCAN for boundary refinement...\n")
dbscan_result <- dbscan::dbscan(tumor_spatial_data[, c("x", "y")], eps = 250, minPts = 30)

# Add the DBSCAN cluster information to the tumor spatial data
tumor_spatial_data$Cluster <- dbscan_result$cluster

# Identify the primary tumor cluster (largest cluster excluding noise)
primary_cluster <- names(which.max(table(tumor_spatial_data$Cluster[tumor_spatial_data$Cluster != 0])))

cat("Primary cluster identified:", primary_cluster, "\n")

# Relabel cells: mark those in the primary cluster as "Tumor" and others as "Non-Tumor"
tumor_spatial_data <- tumor_spatial_data %>%
  mutate(Label = ifelse(Cluster == as.numeric(primary_cluster), "Tumor", "Non-Tumor"))

# Step 5a: Update the CSV with refined labels
cat("Updating CSV with refined tumor labels...\n")
# Update the labels in all_barcodes based on refined tumor boundary
all_barcodes <- all_barcodes %>%
  left_join(
    tumor_spatial_data[, c("cell", "Label")],
    by = c("Barcode" = "cell"),
    suffix = c("_all", "_tumor")
  ) %>%
  mutate(
    Refined_Label = ifelse(!is.na(Label_tumor), Label_tumor, "Non-Tumor")
  ) %>%
  select(Barcode, Refined_Label)

# Save the updated CSV
write_csv(all_barcodes, "JY10_all_barcodes_with_refined_labels.csv")
cat("CSV file 'JY10_all_barcodes_with_refined_labels.csv' generated.\n")

# Step 6: Visualize Refined Tumor Boundary with Outliers Relabeled
cat("Creating refined tumor boundary plot...\n")
refined_boundary_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y, color = Label)) +
  geom_point(size = 0.5) +
  labs(title = "Refined Tumor Boundary with Outliers Relabeled") +
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))

ggsave("refined_tumor_boundary_plot.png", plot = refined_boundary_plot, width = 8, height = 6)
cat("Refined tumor boundary plot saved as 'refined_tumor_boundary_plot.png'.\n")

# Step 7: Load the Annotations File with Cell Types
cat("Loading annotations...\n")
JY10_annotations <- read_csv("JY10_ANNOTATED.csv")

# Ensure the columns are correctly named
if (!all(c("Barcode", "Graph-based") %in% colnames(JY10_annotations))) {
  stop("The columns 'Barcode' and 'Graph-based' must be present in JY10_ANNOTATED.csv")
}

# Standardize barcodes in annotations
JY10_annotations$Barcode <- standardize_barcodes(JY10_annotations$Barcode)

# Step 8: Filter for TILs_T_and_NK_cells
cat("Filtering for TILs_T_and_NK_cells...\n")
target_cells_barcodes <- JY10_annotations %>%
  filter(`Graph-based` == "TILs_T_and_NK_cells") %>%
  select(Barcode)

# Standardize barcodes
target_cells_barcodes$Barcode <- standardize_barcodes(target_cells_barcodes$Barcode)

# Step 9: Merge Target Cells with Spatial Data and Tumor Labels
cat("Merging TILs_T_and_NK_cells with spatial data...\n")
target_spatial_data <- spatial_coords %>%
  filter(cell %in% target_cells_barcodes$Barcode) %>%
  left_join(tumor_spatial_data[, c("cell", "Label")], by = "cell")

# Replace NA in Label with "Non-Tumor" (cells not in the tumor cluster)
target_spatial_data$Label[is.na(target_spatial_data$Label)] <- "Non-Tumor"

# Step 10: Separate Target Cells Inside and Outside Tumor Boundary
cat("Separating TILs_T_and_NK_cells based on tumor boundary...\n")
target_inside_boundary <- target_spatial_data %>%
  filter(Label == "Tumor")

target_outside_boundary <- target_spatial_data %>%
  filter(Label == "Non-Tumor")

# Debugging: Check counts
cat("Number of TILs_T_and_NK_cells inside tumor boundary:", nrow(target_inside_boundary), "\n")
cat("Number of TILs_T_and_NK_cells outside tumor boundary:", nrow(target_outside_boundary), "\n")

# Step 11: Extract Counts Data for Target Cells Only
cat("Extracting counts data for TILs_T_and_NK_cells...\n")
# Get counts data from the Spatial.008um assay
counts_data <- GetAssayData(JY10, assay = "Spatial.008um", slot = "counts")

# Standardize barcodes in counts_data
colnames(counts_data) <- standardize_barcodes(colnames(counts_data))

# Filter counts_data to include only target cells
target_cells <- target_cells_barcodes$Barcode
counts_data_target <- counts_data[, colnames(counts_data) %in% target_cells]

# Confirm dimensions
cat("Counts data dimensions for TILs_T_and_NK_cells:", dim(counts_data_target), "\n")

# Step 12: Identify Inside and Outside Target Cells in Counts Data
cat("Identifying inside and outside TILs_T_and_NK_cells in counts data...\n")
inside_cells <- intersect(target_inside_boundary$cell, colnames(counts_data_target))
outside_cells <- intersect(target_outside_boundary$cell, colnames(counts_data_target))

# Confirm cell counts
cat("Number of inside TILs_T_and_NK_cells:", length(inside_cells), "\n")
cat("Number of outside TILs_T_and_NK_cells:", length(outside_cells), "\n")

# Step 13: Create Seurat Object with Target Cells Only
cat("Creating Seurat object with TILs_T_and_NK_cells only...\n")
seurat_obj <- CreateSeuratObject(counts = counts_data_target)

# Standardize cell names in seurat_obj
seurat_obj <- RenameCells(seurat_obj, new.names = standardize_barcodes(Cells(seurat_obj)))

# Step 14: Assign Group Labels to Seurat Object
cat("Assigning group labels to TILs_T_and_NK_cells...\n")
# Create a named vector for group labels
group_labels <- ifelse(Cells(seurat_obj) %in% inside_cells, "Inside",
                       ifelse(Cells(seurat_obj) %in% outside_cells, "Outside", NA))
names(group_labels) <- Cells(seurat_obj)

# Add the group labels to the metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = group_labels, col.name = "group")

# Verify that 'group' is in metadata
cat("Metadata columns in seurat_obj:\n")
print(colnames(seurat_obj@meta.data))

# Check the distribution of 'group' labels
cat("Distribution of 'group' labels:\n")
print(table(seurat_obj$group, useNA = "ifany"))

# Step 15: Subset Seurat Object to Cells with Group Labels
cat("Subsetting Seurat object to cells with group labels...\n")
# Subset using cell names
cells_to_keep <- rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$group)]
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

# Confirm subsetting
cat("Number of cells after subsetting:", ncol(seurat_obj), "\n")
cat("Group label distribution after subsetting:\n")
print(table(seurat_obj$group))

# Step 16: Proceed with Differential Expression Analysis
cat("Performing differential expression analysis...\n")
# Set the active identity to 'group'
Idents(seurat_obj) <- seurat_obj$group

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# Perform differential expression analysis
de_results <- FindMarkers(seurat_obj, ident.1 = "Inside", ident.2 = "Outside", test.use = "wilcox")

# Debugging: Check the results
cat("Differential expression results (top 10 genes):\n")
print(head(de_results, 10))

# Step 17: Filter for Significant Genes
cat("Filtering for significant genes...\n")
significant_genes <- de_results %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

cat("Number of significant genes:", nrow(significant_genes), "\n")

# Save significant genes to a CSV file
write_csv(significant_genes, "significant_genes_TILs_T_and_NK_cells.csv")
cat("Significant genes saved to 'significant_genes_TILs_T_and_NK_cells.csv'.\n")

# Step 18: Create a Volcano Plot
cat("Creating a volcano plot...\n")

# Prepare data for the volcano plot
volcano_data <- de_results %>%
  mutate(
    gene = rownames(de_results),
    neg_log10_padj = -log10(p_val_adj)
  )

# Basic volcano plot using ggplot2
volcano_plot <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(alpha = 0.6, color = "gray") +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")

# Highlight significant genes
volcano_plot <- volcano_plot +
  geom_point(
    data = volcano_data %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
    aes(x = avg_log2FC, y = neg_log10_padj),
    color = "red",
    alpha = 0.6
  )

# Save the volcano plot
ggsave("volcano_plot_TILs_T_and_NK_cells.png", plot = volcano_plot, width = 8, height = 6)
cat("Volcano plot saved as 'volcano_plot_TILs_T_and_NK_cells.png'.\n")

# Display the plot
print(volcano_plot)

# Step 19: Create a Heatmap of Top Significant Genes
cat("Creating a heatmap of top significant genes...\n")

# Select top 20 significant genes
top_genes <- head(rownames(significant_genes), 20)

# Extract normalized expression data for these genes
expression_data <- GetAssayData(seurat_obj, slot = "data")[top_genes, ]

# Scale the expression data for better visualization
scaled_expression <- t(scale(t(expression_data)))

# Create annotation for columns (cells)
annotation_col <- data.frame(Group = seurat_obj$group)
rownames(annotation_col) <- colnames(expression_data)

# Generate heatmap
heatmap_plot <- pheatmap(
  scaled_expression,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  fontsize_row = 8,
  main = "Heatmap of Top 20 Significant Genes"
)

# Save the heatmap
png("heatmap_top_genes_TILs_T_and_NK_cells.png", width = 800, height = 600)
print(heatmap_plot)
dev.off()
cat("Heatmap saved as 'heatmap_top_genes_TILs_T_and_NK_cells.png'.\n")

# Step 20: Create Violin Plots for Top Significant Genes
cat("Creating violin plots for top significant genes...\n")

# Select top 6 significant genes for visualization
top_violin_genes <- head(rownames(significant_genes), 6)

# Create violin plots
violin_plot <- VlnPlot(
  seurat_obj,
  features = top_violin_genes,
  group.by = "group",
  pt.size = 0.1,
  ncol = 3
) +
  theme_minimal() +
  labs(title = "Violin Plots of Top Significant Genes")

# Save the violin plot
ggsave("violin_plots_top_genes_TILs_T_and_NK_cells.png", plot = violin_plot, width = 12, height = 8)
cat("Violin plots saved as 'violin_plots_top_genes_TILs_T_and_NK_cells.png'.\n")

# Display the plot
print(violin_plot)

cat("Script execution completed with analysis on TILs_T_and_NK_cells.\n")
