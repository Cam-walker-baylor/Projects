# Load Necessary Libraries
library(Seurat)
library(tidyverse)
library(dbscan)
library(ggplot2)
library(sp)        # For spatial operations

# Set options for debugging
options(stringsAsFactors = FALSE)

# Define a function to standardize barcodes
standardize_barcodes <- function(barcodes) {
  barcodes <- trimws(barcodes)    # Remove whitespace
  barcodes <- toupper(barcodes)   # Convert to uppercase
  return(barcodes)
}

# Step 1: Load Data
cat("Loading hypoxia labels...\n")
data <- read_csv("barcodes/JY10_Hypoxia.csv")

# Ensure required columns are present
if (!all(c("Barcode", "Hypoxia") %in% colnames(data))) {
  stop("The columns 'Barcode' and 'Hypoxia' must be present in JY10_Hypoxia.csv")
}

# Standardize barcode formats in data
data$Barcode <- standardize_barcodes(data$Barcode)

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

# Step 3: Merge Hypoxia Labels with Spatial Coordinates
cat("Merging hypoxia labels with spatial coordinates...\n")
tumor_spatial_data <- spatial_coords %>%
  dplyr::left_join(data, by = c("cell" = "Barcode")) %>%
  dplyr::mutate(Label = ifelse(Hypoxia == "Hypoxia", "Hypoxia", "Non-Hypoxia"))

# Debugging: Check number of matched hypoxia cells
cat("Number of hypoxia cells with spatial coordinates:", sum(tumor_spatial_data$Label == "Hypoxia"), "\n")

# Step 4: Initial Tumor Identification Using DBSCAN
cat("Applying DBSCAN for initial tumor identification...\n")
# Extract hypoxia-labeled cells
hypoxia_cells <- tumor_spatial_data %>% dplyr::filter(Label == "Hypoxia")

# Apply DBSCAN clustering
dbscan_result <- dbscan::dbscan(hypoxia_cells[, c("x", "y")], eps = 100, minPts = 5)

# Add cluster information to hypoxia_cells
hypoxia_cells$Cluster <- dbscan_result$cluster

# Identify the primary tumor cluster (largest cluster excluding noise)
primary_cluster <- names(which.max(table(hypoxia_cells$Cluster[hypoxia_cells$Cluster != 0])))
cat("Primary cluster identified:", primary_cluster, "\n")

# Filter hypoxia_cells to only include the primary cluster
primary_tumor_cells <- hypoxia_cells %>%
  dplyr::filter(Cluster == as.numeric(primary_cluster))

# Step 5: Create Convex Hull Around Primary Tumor Cluster
cat("Creating convex hull around the primary tumor cluster...\n")
# Extract coordinates
tumor_coords <- primary_tumor_cells[, c("x", "y")]

# Create a convex hull polygon
hull_indices <- chull(tumor_coords)
hull_coords <- tumor_coords[hull_indices, ]

# Close the polygon by adding the first point at the end
hull_coords <- rbind(hull_coords, hull_coords[1, ])

# Create SpatialPolygon
polygon <- sp::Polygon(hull_coords)
polygons <- sp::Polygons(list(polygon), ID = "1")
tumor_boundary_polygon <- sp::SpatialPolygons(list(polygons))

# Step 6: Label All Cells Within the Polygon as "Tumor"
cat("Relabeling cells within the tumor boundary...\n")
# Create SpatialPoints object for all cells
all_cells_coords <- tumor_spatial_data[, c("x", "y")]
all_cells_points <- sp::SpatialPoints(all_cells_coords)

# Identify cells inside the tumor polygon
inside_tumor <- sp::over(all_cells_points, tumor_boundary_polygon)
tumor_spatial_data$Label <- ifelse(!is.na(inside_tumor), "Tumor", "Non-Tumor")

# Step 7: Visualize the Results
cat("Creating plots...\n")
# Plot initial tumor cells identified by DBSCAN
initial_tumor_plot <- ggplot() +
  geom_point(data = tumor_spatial_data, aes(x = x, y = y), color = "gray", size = 0.5) +
  geom_point(data = primary_tumor_cells, aes(x = x, y = y), color = "red", size = 0.5) +
  labs(title = "Initial Tumor Cells Identified by DBSCAN") +
  theme_minimal()

# Save the initial tumor plot
ggsave("JY10_Tumor/initial_tumor_plot.png", plot = initial_tumor_plot, width = 8, height = 6)
cat("Initial tumor plot saved as 'initial_tumor_plot.png'.\n")

# Plot all cells with updated labels
boundary_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y, color = Label)) +
  geom_point(size = 0.5) +
  labs(title = "Tumor Boundary Using Convex Hull of DBSCAN Cluster") +
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))

# Overlay tumor boundary
boundary_plot <- boundary_plot +
  geom_path(data = as.data.frame(hull_coords), aes(x = x, y = y),
            color = "blue", linetype = "dashed", size = 0.5)

# Save the refined boundary plot
ggsave("JY10_Tumor/tumor_boundary_plot.png", plot = boundary_plot, width = 8, height = 6)
cat("Tumor boundary plot saved as 'tumor_boundary_plot.png'.\n")

# Step 8: Save the Relabeled Barcodes
cat("Saving relabeled barcodes...\n")
# Select the cell barcodes and their updated labels
relabeled_barcodes <- tumor_spatial_data %>% dplyr::select(cell, Label)

# Save to a CSV file
write_csv(relabeled_barcodes, "JY10_Tumor/relabeled_barcodes.csv")
cat("Relabeled barcodes saved as 'relabeled_barcodes.csv'.\n")

# Optional: Display the plots in your R session
print(initial_tumor_plot)
print(boundary_plot)
