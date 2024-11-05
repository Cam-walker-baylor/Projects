# Load necessary libraries
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(dbscan)

# Step 1: Load Data
# Load the data file containing barcodes labeled as "Tumor"
data <- read_csv("JY10_Tumor_Micro_ENV.csv")

# Filter for barcodes labeled as "Tumor"
tumor_barcodes <- data %>%
  filter(`2` == "Tumor") %>%
  select(Barcode)

# Step 2: Load Seurat Object and Extract Spatial Coordinates
# Load the Seurat object containing spatial data
JY10 <- readRDS("rds_files/JY10_object.rds")

# Extract spatial coordinates for the "slice1.008um" image
spatial_coords <- GetTissueCoordinates(JY10, image = "slice1.008um")

# Step 3: Merge Tumor Labels with Spatial Coordinates
# Merge the filtered tumor barcodes with spatial coordinates
tumor_spatial_data <- spatial_coords %>%
  filter(cell %in% tumor_barcodes$Barcode) %>%
  mutate(Label = "Tumor")

# Step 4: Visualize Initial Tumor-Labeled Cells
# Plot the initial distribution of tumor-labeled cells
initial_tumor_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y)) +
  geom_point(color = "red", size = 0.5) +
  labs(title = "Spatial Distribution of Tumor-Labeled Cells") +
  theme_minimal()

# Save the initial tumor plot to a file
ggsave("initial_tumor_plot.png", plot = initial_tumor_plot, width = 8, height = 6)

# Step 5: Apply DBSCAN for Boundary Refinement
# Run DBSCAN on the tumor-labeled spatial coordinates
# Adjust eps and minPts to refine the boundary further
dbscan_result <- dbscan::dbscan(tumor_spatial_data[, c("x", "y")], eps = 250, minPts = 30)

# Add the DBSCAN cluster information to the tumor spatial data
tumor_spatial_data$Cluster <- dbscan_result$cluster

# Identify the primary tumor cluster (largest cluster with most points)
primary_cluster <- names(which.max(table(tumor_spatial_data$Cluster[tumor_spatial_data$Cluster != 0])))

# Relabel cells: mark those in the primary cluster as "Tumor" and outliers as "Non-Tumor"
tumor_spatial_data <- tumor_spatial_data %>%
  mutate(Label = ifelse(Cluster == as.numeric(primary_cluster), "Tumor", "Non-Tumor"))

# Step 6: Visualize Refined Tumor Boundary with Outliers Relabeled
# Plot the refined tumor boundary, marking outliers in gray
refined_boundary_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y, color = Label)) +
  geom_point(size = 0.5) +
  labs(title = "More Refined Tumor Boundary with Outliers Relabeled") +
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))

# Save the refined boundary plot to a file
ggsave("more_refined_tumor_boundary_plot.png", plot = refined_boundary_plot, width = 8, height = 6)

# Step 7: Load Myeloid Annotations
# Load the annotations file with cell types
JY10_annotations <- read_csv("JY10_ANNOTATED.csv")

# Filter for Myeloid cells
myeloid_barcodes <- JY10_annotations %>%
  filter(`Graph-based` == "Myeloid_cells") %>%
  select(Barcode)

# Merge with spatial data to add tumor boundary labels
# Join myeloid barcodes with spatial data
myeloid_spatial_data <- spatial_coords %>%
  filter(cell %in% myeloid_barcodes$Barcode) %>%
  left_join(tumor_spatial_data, by = c("cell" = "cell"))

# Separate Myeloid cells into inside and outside tumor boundary
# "Inside" if Label is "Tumor", "Outside" if Label is "Non-Tumor"
myeloid_inside_boundary <- myeloid_spatial_data %>%
  filter(Label == "Tumor")

myeloid_outside_boundary <- myeloid_spatial_data %>%
  filter(Label == "Non-Tumor")

# Step 8: Access Counts Data for Differential Expression
# Access the counts layer directly from the Spatial.008um assay
counts_data <- JY10[["Spatial.008um"]]@layers[["counts"]]

# Set row and column names if not already present
rownames(counts_data) <- rownames(JY10)  # Replace with gene names if available
colnames(counts_data) <- colnames(JY10)  # Replace with cell names if available

# Step 9: Perform Differential Expression Analysis
# Split Cells by Group
inside_cells <- colnames(counts_data)[JY10_myeloid$boundary_status == "Inside"]
outside_cells <- colnames(counts_data)[JY10_myeloid$boundary_status == "Outside"]

# Compute Average Expression and Perform Statistical Tests
de_results <- data.frame(
  gene = rownames(counts_data),
  avg_inside = rowMeans(counts_data[, inside_cells, drop = FALSE]),
  avg_outside = rowMeans(counts_data[, outside_cells, drop = FALSE]),
  p_value = apply(counts_data, 1, function(x) {
    wilcox.test(x[inside_cells], x[outside_cells])$p.value
  })
)

# Calculate log2 fold change
de_results <- de_results %>%
  mutate(
    log2FC = log2(avg_inside + 1) - log2(avg_outside + 1),
    p_adj = p.adjust(p_value, method = "fdr")
  )

# Filter for significant genes (adjust threshold as needed)
significant_genes <- de_results %>%
  filter(p_adj < 0.05 & abs(log2FC) > 0.5)

# Step 10: Visualization

# Volcano Plot
volcano_plot <- ggplot(de_results, aes(x = log2FC, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05 & abs(log2FC) > 0.5)) +
  labs(
    title = "Differential Expression: Inside vs. Outside Tumor Boundary (Myeloid Cells)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"))

# Save volcano plot
ggsave("volcano_plot_myeloid_inside_vs_outside.png", plot = volcano_plot, width = 8, height = 6)

# Heatmap for Top Genes
top_genes <- significant_genes %>%
  arrange(-abs(log2FC)) %>%
  head(20) %>%
  pull(gene)

# Prepare data for heatmap
heatmap_data <- counts_data[top_genes, c(inside_cells, outside_cells)]

# Standardize gene expression for visualization
heatmap_data <- t(scale(t(as.matrix(heatmap_data))))

# Convert to long format for ggplot
heatmap_long <- as.data.frame(heatmap_data) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
  mutate(group = ifelse(cell %in% inside_cells, "Inside", "Outside"))

# Heatmap plot
heatmap_plot <- ggplot(heatmap_long, aes(x = cell, y = gene, fill = expression)) +
  geom_tile() +
  facet_wrap(~group, scales = "free_x", ncol = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(
    title = "Top Differentially Expressed Genes (Myeloid Cells)",
    x = "Cells",
    y = "Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Save heatmap plot
ggsave("heatmap_myeloid_inside_vs_outside.png", plot = heatmap_plot, width = 10, height = 8)