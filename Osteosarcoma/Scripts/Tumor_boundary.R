# Load necessary libraries
library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(dbscan)

# Step 1: Load Data
# Load the data file containing barcodes labeled as "Tumor"
data <- read_csv("barcodes/JY10_Tumor_Micro_ENV.csv")

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
ggsave("JY10_Tumor/initial_tumor_plot.png", plot = initial_tumor_plot, width = 8, height = 6)

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
ggsave("JY10_Tumor/more_refined_tumor_boundary_plot.png", plot = refined_boundary_plot, width = 8, height = 6)

# Display final refined boundary plot
refined_boundary_plot

# Load the annotations file with cell types
JY10_annotations <- read_csv("JY10_ANNOTATED.csv")

# Step 1: Filter for Myeloid cells
myeloid_barcodes <- JY10_annotations %>%
  filter(`Graph-based` == "Myeloid_cells") %>%
  select(Barcode)

# Step 2: Merge with spatial data to add tumor boundary labels
# Join myeloid barcodes with spatial data
myeloid_spatial_data <- spatial_coords %>%
  filter(cell %in% myeloid_barcodes$Barcode) %>%
  left_join(tumor_spatial_data, by = c("cell" = "cell"))

# Step 3: Separate Myeloid cells into inside and outside tumor boundary
# "Inside" if Label is "Tumor", "Outside" if Label is "Non-Tumor"
myeloid_inside_boundary <- myeloid_spatial_data %>%
  filter(Label == "Tumor")

myeloid_outside_boundary <- myeloid_spatial_data %>%
  filter(Label == "Non-Tumor")

# Output the results
# You can view or save these groups as needed
print("Myeloid Cells Inside Tumor Boundary:")
print(myeloid_inside_boundary)

print("Myeloid Cells Outside Tumor Boundary:")
print(myeloid_outside_boundary)

# Optional: Save to CSV for further use
write_csv(myeloid_inside_boundary, "myeloid_inside_boundary.csv")
write_csv(myeloid_outside_boundary, "myeloid_outside_boundary.csv")

####################### NEW

# Extract counts data from the Spatial.008um assay
counts_data <- JY10[["Spatial.008um"]]@layers[["counts"]]

# Assign cell names to counts_data if they are missing
if (is.null(colnames(counts_data))) {
  colnames(counts_data) <- Cells(JY10[["Spatial.008um"]])
}

# Verify that counts_data now has cell names
cat("Column names in counts_data after assignment:\n")
print(head(colnames(counts_data)))

# Identify cells that match counts_data for Inside and Outside groups
inside_cells <- intersect(myeloid_inside_boundary$cell, colnames(counts_data))
outside_cells <- intersect(myeloid_outside_boundary$cell, colnames(counts_data))

# Confirm cell counts
cat("Number of Inside Cells:", length(inside_cells), "\n")
cat("Number of Outside Cells:", length(outside_cells), "\n")

# Subset counts_data for Inside and Outside Cells
counts_inside <- counts_data[, inside_cells, drop = FALSE]
counts_outside <- counts_data[, outside_cells, drop = FALSE]

# Perform differential expression analysis
de_results <- data.frame(
  gene = rownames(counts_data),
  avg_inside = rowMeans(counts_inside),
  avg_outside = rowMeans(counts_outside),
  p_value = apply(counts_data, 1, function(x) {
    wilcox.test(x[inside_cells], x[outside_cells])$p.value
  })
)

# Calculate log2 fold change and adjust p-values
de_results <- de_results %>%
  mutate(
    log2FC = log2(avg_inside + 1) - log2(avg_outside + 1),
    p_adj = p.adjust(p_value, method = "fdr")
  )

# Filter for significant genes
significant_genes <- de_results %>%
  filter(p_adj < 0.05 & abs(log2FC) > 0.5)






################ INspect Data

# View the overall structure of the Seurat object
str(JY10)

# List available assays
Assays(JY10)

# View structure of the Spatial.008um assay
str(JY10[["Spatial.008um"]])

# Check available layers within the Spatial.008um assay
Layers(JY10[["Spatial.008um"]])

# Extract and examine the counts layer
counts_data <- JY10[["Spatial.008um"]]@layers[["counts"]]

# Check dimensions
dim(counts_data)

# View row names (genes) and column names (cells) if available
head(rownames(counts_data))
head(colnames(counts_data))

# View the first few rows of metadata
head(JY10@meta.data)

# Check structure of metadata
str(JY10@meta.data)

# List available images in JY10
Images(JY10)

# Extract spatial coordinates from a specific image, e.g., "slice1.008um"
spatial_coords <- GetTissueCoordinates(JY10, image = "slice1.008um")
head(spatial_coords)

# List reductions available in JY10
Reductions(JY10)


