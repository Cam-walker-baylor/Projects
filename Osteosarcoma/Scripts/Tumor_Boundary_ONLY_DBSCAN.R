# Load Necessary Libraries
library(Seurat)         # For handling Seurat objects
library(tidyverse)      # Includes readr, dplyr, ggplot2, etc.
library(dbscan)         # For density-based clustering
library(ggplot2)        # For plotting

# Set options for debugging
options(stringsAsFactors = FALSE)

# Define a function to standardize barcodes
standardize_barcodes <- function(barcodes) {
  barcodes <- trimws(barcodes)    # Remove whitespace
  barcodes <- toupper(barcodes)   # Convert to uppercase
  return(barcodes)
}

# Step 1: Load Data
cat("Loading tumor microenvironment labels...\n")
data <- read_csv("barcodes/JY10_Tumor_Micro_ENV.csv")

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

# Step 4: Visualize Initial Tumor-Labeled Cells
cat("Creating initial tumor plot...\n")
initial_tumor_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y)) +
  geom_point(color = "red", size = 0.5) +
  labs(title = "Spatial Distribution of Tumor-Labeled Cells") +
  theme_minimal()

# Save the initial tumor plot to a file
ggsave("JY10_Tumor/initial_tumor_plot.png", plot = initial_tumor_plot, width = 8, height = 6)
cat("Initial tumor plot saved as 'initial_tumor_plot.png'.\n")

# Step 5: Apply DBSCAN for Boundary Refinement
cat("Applying DBSCAN for boundary refinement...\n")
# Adjust 'eps' and 'minPts' as needed for your data
dbscan_result <- dbscan::dbscan(tumor_spatial_data[, c("x", "y")], eps = 250, minPts = 15)

# Add the DBSCAN cluster information to the tumor spatial data
tumor_spatial_data$Cluster <- dbscan_result$cluster

# Identify the primary tumor cluster (largest cluster excluding noise)
primary_cluster <- names(which.max(table(tumor_spatial_data$Cluster[tumor_spatial_data$Cluster != 0])))

cat("Primary cluster identified:", primary_cluster, "\n")

# Relabel cells: mark those in the primary cluster as "Tumor" and others as "Non-Tumor"
tumor_spatial_data <- tumor_spatial_data %>%
  mutate(Label = ifelse(Cluster == as.numeric(primary_cluster), "Tumor", "Non-Tumor"))

# Step 6: Visualize Refined Tumor Boundary with Outliers Relabeled
cat("Creating refined tumor boundary plot...\n")
refined_boundary_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y, color = Label)) +
  geom_point(size = 0.5) +
  labs(title = "Refined Tumor Boundary with Outliers Relabeled") +
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))

# Save the refined boundary plot to a file
ggsave("JY10_Tumor/refined_tumor_boundary_plot.png", plot = refined_boundary_plot, width = 8, height = 6)
cat("Refined tumor boundary plot saved as 'refined_tumor_boundary_plot.png'.\n")

# Step 7: Display the Plots
# (Optional) Display the plots in your R session
print(initial_tumor_plot)
print(refined_boundary_plot)
