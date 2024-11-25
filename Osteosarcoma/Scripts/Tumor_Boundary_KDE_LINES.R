# Load necessary libraries
library(readr)     # For reading CSV files
library(dplyr)     # For data manipulation
library(Seurat)    # For working with Seurat objects
library(ggplot2)   # For plotting
library(MASS)      # For KDE using kde2d
library(sp)        # For spatial operations
library(ggpubr)    # For plotting
library(dbscan)    # For HDBSCAN clustering

# Step 1: Load and Prepare the Data

# Load hypoxia-labeled barcodes data
data <- read_csv("barcodes/JY10_Hypoxia2.csv")

# Load the Seurat object containing spatial data
JY10 <- readRDS("rds_files/JY10_object.rds")

# Extract spatial coordinates for the "slice1.008um" image
spatial_coords <- GetTissueCoordinates(JY10, image = "slice1.008um")

# Merge hypoxia labels with spatial coordinates
tumor_spatial_data <- spatial_coords %>%
  left_join(data %>% dplyr::select(Barcode, Hypoxia), by = c("cell" = "Barcode")) %>%
  mutate(Label = ifelse(!is.na(Hypoxia), "Hypoxia", "Non-Hypoxia"))

# Step 2: Perform KDE with Optimized Bandwidth and Increased Resolution

# Extract coordinates of hypoxia-labeled cells
hypoxia_cells <- tumor_spatial_data %>% filter(Label == "Hypoxia")
hypoxia_coords <- hypoxia_cells[, c("x", "y")]

# Estimate optimal bandwidths
bandwidth_x <- MASS::bandwidth.nrd(hypoxia_coords$x) / 2  # Adjust divisor as needed
bandwidth_y <- MASS::bandwidth.nrd(hypoxia_coords$y) / 2

# Perform KDE with specified bandwidths and higher resolution
kde_result <- MASS::kde2d(
  hypoxia_coords$x,
  hypoxia_coords$y,
  h = c(bandwidth_x, bandwidth_y),
  n = 1000  # Increased grid resolution
)

# Step 3: Visualize Density Distribution

# Flatten density matrix to a vector
density_values <- as.vector(kde_result$z)

# Plot density distribution
density_df <- data.frame(Density = density_values)
ggplot(density_df, aes(x = Density)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  labs(title = "Density Distribution of Hypoxia-Labeled Cells",
       x = "Density", y = "Frequency") +
  theme_minimal()

# Step 4: Extract Multiple Contours at Different Density Levels

# Define a range of density thresholds (e.g., from 30th to 70th percentile)
density_levels <- quantile(density_values, probs = seq(0.5, 0.9, by = 0.1))

# Initialize list to store contours
contour_list <- list()

# Extract contours at each density level
for (level in density_levels) {
  contours <- contourLines(kde_result$x, kde_result$y, kde_result$z, levels = level)
  contour_list <- c(contour_list, contours)
}

# Step 5: Create Polygons from Contours and Filter by Area

# Create polygons for each contour
tumor_polygons_list <- list()
polygon_id <- 1

for (contour in contour_list) {
  coords <- data.frame(x = contour$x, y = contour$y)
  # Close the polygon if necessary
  if (!all(coords[1, ] == coords[nrow(coords), ])) {
    coords <- rbind(coords, coords[1, ])
  }
  # Create polygon and assign ID
  polygon <- sp::Polygon(coords)
  polygons <- sp::Polygons(list(polygon), ID = as.character(polygon_id))
  tumor_polygons_list[[polygon_id]] <- polygons
  polygon_id <- polygon_id + 1
}

# Combine polygons into SpatialPolygons object
tumor_spatial_polygons <- sp::SpatialPolygons(tumor_polygons_list)

# Filter polygons by area to remove small, insignificant regions
min_area <- 500  # Adjust based on data scale
tumor_polygons_list_filtered <- Filter(function(polygon) {
  area <- polygon@Polygons[[1]]@area
  area >= min_area
}, tumor_polygons_list)

# Recreate SpatialPolygons object with filtered polygons
tumor_spatial_polygons <- sp::SpatialPolygons(tumor_polygons_list_filtered)

# Step 6: Label All Cells Within the Polygons as "Tumor"

# Create SpatialPoints object for all cells
all_cells_coords <- tumor_spatial_data[, c("x", "y")]
all_cells_points <- sp::SpatialPoints(all_cells_coords)

# Identify cells inside any of the tumor polygons
inside_tumor <- sp::over(all_cells_points, tumor_spatial_polygons)
tumor_spatial_data$Label <- ifelse(!is.na(inside_tumor), "Tumor", "Non-Tumor")

# Step 7: Visualize the Results

# Plot all cells with updated labels
boundary_plot <- ggplot(tumor_spatial_data, aes(x = x, y = y, color = Label)) +
  geom_point(size = 0.5) +
  labs(title = "Improved Tumor Boundary Based on Multiple Contours") +
  theme_minimal() +
  scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))

# Overlay tumor boundaries
for (polygon in tumor_polygons_list_filtered) {
  coords <- polygon@Polygons[[1]]@coords
  boundary_plot <- boundary_plot +
    geom_polygon(data = data.frame(x = coords[, 1], y = coords[, 2]),
                 aes(x = x, y = y),
                 fill = NA, color = "blue", linetype = "dashed", size = 0.5)
}

# Save the plot
ggsave("JY10_Tumor/improved_tumor_boundary.png", plot = boundary_plot, width = 8, height = 6)

# Step 8: Save the Relabeled Barcodes

# Select the cell barcodes and their updated labels
relabeled_barcodes <- tumor_spatial_data %>% dplyr::select(cell, Label)

# Save to a CSV file
write_csv(relabeled_barcodes, "relabeled_barcodes.csv")
