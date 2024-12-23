library(dplyr)
library(ggplot2)
library(Seurat)
library(Matrix)

#' Configuration settings for binning analysis
#' @param bin_size Bin size (default 11 for 55µm)
#' @param assay Assay name
#' @return Named list of parameters
get_config <- function(bin_size = 11, assay = "Spatial.008um") {
  list(
    bin_size = bin_size,
    assay = assay,
    required_packages = c("dplyr", "ggplot2", "Seurat", "Matrix")
  )
}

#' Create directory structure
#' @param base_dir Base directory
#' @return Named list of paths
create_dirs <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    plots = file.path(base_dir, "plots"),
    data = file.path(base_dir, "data"),
    stats = file.path(base_dir, "stats")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

#' Create spatial bins
#' @param coords Coordinates data frame
#' @param counts Count matrix
#' @param params Analysis parameters
#' @return List of binned data
bin_spatial_data <- function(coords, counts, params) {
  bin_factor <- params$bin_size * params$scale_factor
  bin_x <- floor(coords$x / bin_factor)
  bin_y <- floor(coords$y / bin_factor)
  bin_id <- paste0("bin_", bin_x, "_", bin_y)
  
  binned_counts <- sapply(split(seq_len(ncol(counts)), bin_id), function(idx) {
    if(length(idx) > 0) Matrix::rowSums(counts[, idx, drop = FALSE]) else rep(0, nrow(counts))
  })
  
  bin_coords <- data.frame(
    bin_id = unique(bin_id),
    x = tapply(coords$x, bin_id, mean),
    y = tapply(coords$y, bin_id, mean)
  )
  
  list(counts = binned_counts, coords = bin_coords, bin_ids = bin_id)
}

#' Generate QC metrics
#' @param original Original count matrix
#' @param binned Binned data
#' @param sample_id Sample identifier
#' @return Data frame of statistics
generate_qc_metrics <- function(original, binned, sample_id) {
  total_counts <- colSums(binned$counts)
  data.frame(
    sample_id = sample_id,
    original_spots = ncol(original),
    binned_spots = ncol(binned$counts),
    median_counts = median(total_counts),
    mean_counts = mean(total_counts),
    genes = nrow(binned$counts),
    timestamp = Sys.time()
  )
}

#' Create visualization
#' @param object Seurat object
#' @param title Plot title
#' @return ggplot object
visualize_bins <- function(object, title) {
  ggplot(object@meta.data, aes(x = x, y = y, color = nCount_RNA)) +
    geom_point(size = 2) +
    scale_color_viridis_c() +
    theme_minimal() +
    ggtitle(title) +
    theme(aspect.ratio = 1)
}

#' Process single sample
#' @param object Seurat object
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param config Configuration parameters
#' @return List of results
process_sample <- function(object, sample_id, output_dir, config) {
  dirs <- create_dirs(file.path(output_dir, sample_id))
  
  DefaultAssay(object) <- config$assay
  coords <- GetTissueCoordinates(object)
  counts <- GetAssayData(object, slot = "counts")
  config$scale_factor <- object@images[[paste0("slice1.", sub("Spatial.", "", config$assay))]]@scale.factors$spot
  
  binned <- bin_spatial_data(coords, counts, config)
  binned_obj <- CreateSeuratObject(counts = binned$counts)
  binned_obj@meta.data[c("x", "y")] <- binned$coords[match(colnames(binned_obj), binned$coords$bin_id), c("x", "y")]
  
  plot <- visualize_bins(binned_obj, paste("Total Counts per Bin -", sample_id))
  stats <- generate_qc_metrics(counts, binned, sample_id)
  
  ggsave(file.path(dirs$plots, "spatial_distribution.pdf"), plot, width = 10, height = 10)
  write.csv(stats, file.path(dirs$stats, "statistics.csv"), row.names = FALSE)
  saveRDS(binned_obj, file.path(dirs$data, "binned_object.rds"))
  
  list(object = binned_obj, stats = stats, plot = plot)
}

#' Main pipeline function
#' @param samples Named list of Seurat objects
#' @param output_dir Output directory
#' @param bin_size Bin size
#' @param assay Assay name
#' @return Named list of results
run_binning_pipeline <- function(samples, output_dir, bin_size = 11, assay = "Spatial.008um") {
  config <- get_config(bin_size, assay)
  results <- lapply(names(samples), function(sample_id) {
    process_sample(samples[[sample_id]], sample_id, output_dir, config)
  })
  names(results) <- names(samples)
  results
}

# Example usage:
samples <- list(
  JY4 = readRDS("rds_files/JY4_object.rds"),
  JY8 = readRDS("rds_files/JY8_object.rds"),
  JY10 = readRDS("rds_files/JY10_object.rds"),
  JY23 = readRDS("rds_files/JY23_object.rds")
)

results <- run_binning_pipeline(samples, "visium_analysis")