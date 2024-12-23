library(readr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(MASS)
library(sp)
library(ggpubr)
library(dbscan)

#' Prepare spatial data with hypoxia labels
#' @param hypoxia_file Path to hypoxia CSV
#' @param seurat_object Seurat object
#' @return Combined spatial data frame
prepare_spatial_data <- function(hypoxia_file, seurat_object) {
  hypoxia_data <- read_csv(hypoxia_file)
  spatial_coords <- GetTissueCoordinates(seurat_object, image = "slice1.008um")
  
  spatial_coords %>%
    left_join(hypoxia_data %>% select(Barcode, Hypoxia), 
              by = c("cell" = "Barcode")) %>%
    mutate(Label = ifelse(!is.na(Hypoxia), "Hypoxia", "Non-Hypoxia"))
}

#' Perform KDE analysis
#' @param coords Coordinates data frame
#' @param bandwidth_divisor Divisor for bandwidth calculation
#' @param grid_size Grid resolution
#' @return KDE results list
perform_kde <- function(coords, bandwidth_divisor = 2, grid_size = 500) {
  bandwidth_x <- MASS::bandwidth.nrd(coords$x) / bandwidth_divisor
  bandwidth_y <- MASS::bandwidth.nrd(coords$y) / bandwidth_divisor
  
  MASS::kde2d(
    coords$x,
    coords$y,
    h = c(bandwidth_x, bandwidth_y),
    n = grid_size
  )
}

#' Extract contours from KDE results
#' @param kde_result KDE results
#' @param probs Probability range for density thresholds
#' @return List of contours
extract_contours <- function(kde_result, probs = seq(0.3, 0.7, by = 0.1)) {
  density_values <- as.vector(kde_result$z)
  density_levels <- quantile(density_values, probs = probs)
  
  contour_list <- list()
  for (level in density_levels) {
    contours <- contourLines(kde_result$x, kde_result$y, kde_result$z, 
                             levels = level)
    contour_list <- c(contour_list, contours)
  }
  contour_list
}

#' Create and filter polygons
#' @param contour_list List of contours
#' @param min_area Minimum area threshold
#' @return Filtered SpatialPolygons object
create_filtered_polygons <- function(contour_list, min_area = 500) {
  polygon_list <- lapply(seq_along(contour_list), function(i) {
    contour <- contour_list[[i]]
    coords <- data.frame(x = contour$x, y = contour$y)
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
    sp::Polygons(list(sp::Polygon(coords)), ID = as.character(i))
  })
  
  filtered_polygons <- Filter(function(polygon) {
    polygon@Polygons[[1]]@area >= min_area
  }, polygon_list)
  
  sp::SpatialPolygons(filtered_polygons)
}

#' Label cells based on polygon boundaries
#' @param spatial_data Spatial data frame
#' @param polygons SpatialPolygons object
#' @return Updated spatial data frame
label_cells <- function(spatial_data, polygons) {
  points <- sp::SpatialPoints(spatial_data[, c("x", "y")])
  inside <- sp::over(points, polygons)
  spatial_data %>%
    mutate(Label = ifelse(!is.na(inside), "Tumor", "Non-Tumor"))
}

#' Create boundary plot
#' @param spatial_data Labeled spatial data
#' @param polygons Tumor boundary polygons
#' @return ggplot object
create_boundary_plot <- function(spatial_data, polygons) {
  plot <- ggplot(spatial_data, aes(x = x, y = y, color = Label)) +
    geom_point(size = 0.5) +
    labs(title = "Tumor Boundary Detection") +
    theme_minimal() +
    scale_color_manual(values = c("Tumor" = "red", "Non-Tumor" = "gray"))
  
  for (polygon in polygons@polygons) {
    coords <- polygon@Polygons[[1]]@coords
    plot <- plot +
      geom_polygon(data = data.frame(x = coords[, 1], y = coords[, 2]),
                   aes(x = x, y = y),
                   fill = NA, color = "blue", linetype = "dashed", size = 0.5)
  }
  plot
}

#' Main pipeline function
#' @param hypoxia_file Path to hypoxia CSV
#' @param seurat_file Path to Seurat object RDS
#' @param output_prefix Prefix for output files
run_tumor_boundary_pipeline <- function(hypoxia_file, seurat_file, output_prefix) {
  # Load data
  seurat_object <- readRDS(seurat_file)
  spatial_data <- prepare_spatial_data(hypoxia_file, seurat_object)
  
  # Process hypoxia cells
  hypoxia_cells <- spatial_data %>% filter(Label == "Hypoxia")
  kde_result <- perform_kde(hypoxia_cells[, c("x", "y")])
  
  # Generate and process contours
  contours <- extract_contours(kde_result)
  polygons <- create_filtered_polygons(contours)
  
  # Label cells and create visualization
  labeled_data <- label_cells(spatial_data, polygons)
  boundary_plot <- create_boundary_plot(labeled_data, polygons)
  
  # Save outputs
  ggsave(paste0(output_prefix, "_boundary.png"), 
         plot = boundary_plot, width = 8, height = 6)
  write_csv(labeled_data %>% select(cell, Label),
            paste0(output_prefix, "_labels.csv"))
  
  list(
    data = labeled_data,
    plot = boundary_plot,
    polygons = polygons
  )
}

# Usage example:
# results <- run_tumor_boundary_pipeline(
#   "JY23_Hypoxia.csv",
#   "rds_files/JY23_object.rds",
#   "JY23"
# )