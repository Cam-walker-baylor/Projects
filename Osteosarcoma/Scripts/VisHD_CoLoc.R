# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(igraph)
library(reshape2)

# Clear the environment to avoid variable conflicts
rm(list = ls())

# ----------------------------------------------
# **Step 1: Set Parameters and Create Directories**
# ----------------------------------------------

# Set the sample ID
# Set sample ID
sample_id <- "JY4"  # Change this to "JY8", "JY10", or "JY23" as needed

# Set analysis parameters
panel_size <- 100     # Adjust this value as needed
overlap <- 0.5         # Overlap between panels (e.g., 50% overlap)
p_value_threshold <- 0.001  # P-value threshold for significance

# Convert p-value to desired format (removing decimal and leading zeros)
p_value_str <- gsub("^0\\.", "", sprintf("%.3f", p_value_threshold))

# Define output directories with the new p-value format
output_dir <- paste0("outs/", sample_id, "_p", panel_size, "_o", overlap * 100,
                     "_pval_", p_value_str, "/")
plot_dir <- paste0(output_dir, "plots/")
spatial_plot_dir <- paste0(plot_dir, "spatial_plots/")
file_dir <- paste0(output_dir, "files/")

# Create output directories if they don't exist
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

# Print the created directory path to verify
print(output_dir)

# ----------------------------------------------
# **Step 2: Load Data and Prepare Seurat Object**
# ----------------------------------------------

# Load the Seurat object for the sample
seurat_object <- readRDS(paste0("rds_files/", sample_id, "_object.rds"))

# Get the barcodes (cell names) from the Seurat object
barcodes <- Cells(seurat_object)

# Read in the cluster assignments CSV file for the sample
cluster_assignments <- read.csv(paste0(sample_id, "_ANNOTATED.csv"), header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------
# **Step 3: Ensure Barcodes Match Exactly**
# ------------------------------------------------

# Check if barcodes in the CSV are present in the Seurat object
matching_barcodes <- cluster_assignments$Barcode %in% barcodes

# Number of matching barcodes
num_matching_barcodes <- sum(matching_barcodes)
print(paste("Number of matching barcodes:", num_matching_barcodes))

# Keep only the cluster assignments that match the Seurat object barcodes
cluster_assignments <- cluster_assignments[matching_barcodes, ]

# Set barcodes as row names
rownames(cluster_assignments) <- cluster_assignments$Barcode

# Remove the Barcode column to prevent duplication
cluster_assignments$Barcode <- NULL

# ---------------------------------------------------------------
# **Step 4: Add Cluster Assignments to Seurat Object's Metadata**
# ---------------------------------------------------------------

# Ensure the order of barcodes matches between metadata and Seurat object
cluster_assignments <- cluster_assignments[barcodes, , drop = FALSE]

# Add the cluster assignments to the Seurat object's metadata
seurat_object <- AddMetaData(object = seurat_object, metadata = cluster_assignments)

# Assign "Unassigned" to NA cluster assignments
seurat_object@meta.data$Graph.based[is.na(seurat_object@meta.data$Graph.based)] <- "Unassigned"

# ---------------------------------------------
# **Step 5: Extract Spatial Coordinates and Cluster Information**
# ---------------------------------------------

# Specify the image name (adjust if necessary)
image_name <- "slice1.008um"  # Adjust if needed per sample

# Get tissue coordinates for the specified image
coords <- GetTissueCoordinates(seurat_object, image = image_name)

# Add cluster information to the coordinates data frame
coords$cluster <- seurat_object@meta.data[rownames(coords), "Graph.based"]

# Convert cluster to factor for consistency
coords$cluster <- as.factor(coords$cluster)

# **Exclude 'Unassigned' cells from analysis, but keep for spatial positioning**

# Create a subset of coords excluding 'Unassigned' for analysis
coords_for_analysis <- coords %>%
  filter(cluster != "Unassigned")

# -----------------------------------------------------
# **Step 6: Define Grid Parameters for Spatial Analysis**
# -----------------------------------------------------

# Calculate stride based on overlap
stride <- panel_size * (1 - overlap)

# -------------------------------------------
# **Step 7: Generate Grid Coordinates**
# -------------------------------------------

# Determine the range of x and y coordinates
x_min <- min(coords$x)
x_max <- max(coords$x)
y_min <- min(coords$y)
y_max <- max(coords$y)

# Create sequences of x and y positions for the grid
x_seq <- seq(x_min, x_max - panel_size, by = stride)
y_seq <- seq(y_min, y_max - panel_size, by = stride)

# ----------------------------------------------
# **Step 8: Iterate Over Each Grid Panel**
# ----------------------------------------------

# Initialize a list to store results
results_list <- list()
panel_id <- 1

for (x_start in x_seq) {
  for (y_start in y_seq) {
    x_end <- x_start + panel_size
    y_end <- y_start + panel_size
    
    # Select cells within the current panel (include 'Unassigned' only for spatial positioning)
    in_panel <- coords %>%
      filter(x >= x_start & x < x_end & y >= y_start & y < y_end)
    
    # Proceed only if there are cells in the panel
    if (nrow(in_panel) > 0) {
      # Exclude 'Unassigned' from analysis
      in_panel_analysis <- in_panel %>%
        filter(cluster != "Unassigned")
      
      # Proceed only if there are assigned cells in the panel
      if (nrow(in_panel_analysis) > 0) {
        # Calculate the proportion of each cluster
        cluster_counts <- in_panel_analysis %>%
          group_by(cluster) %>%
          summarise(count = n())
        
        total_counts <- sum(cluster_counts$count)
        cluster_counts <- cluster_counts %>%
          mutate(proportion = count / total_counts)
        
        # Add panel information
        cluster_counts <- cluster_counts %>%
          mutate(panel_id = panel_id,
                 x_start = x_start,
                 x_end = x_end,
                 y_start = y_start,
                 y_end = y_end)
        
        # Store the results
        results_list[[panel_id]] <- cluster_counts
        
        # Increment panel ID
        panel_id <- panel_id + 1
      }
    }
  }
}

# Combine the results into a data frame
results_df <- bind_rows(results_list)

# Reset row names
rownames(results_df) <- NULL

# --------------------------------------------
# **Step 9: Visualization of Cluster Proportions**
# --------------------------------------------

# Example: Plot the proportion of a cluster across panels
cluster_of_interest <- "Myeloid_cells"  # Adjust this to your cluster of interest

# Filter results for the cluster of interest
plot_data <- results_df %>%
  filter(cluster == cluster_of_interest)

# Check if plot_data is not empty
if (nrow(plot_data) > 0) {
  # Create a plot with panel positions
  my_plot <- ggplot(plot_data, aes(x = x_start, y = y_start, fill = proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    labs(title = paste("Proportion of", cluster_of_interest, "in", sample_id,
                       "(Panel:", panel_size, "um,", overlap * 100, "% overlap)"),
         x = "X Coordinate",
         y = "Y Coordinate",
         fill = "Proportion") +
    theme_minimal()
  
  # Save the plot
  plot_filename <- paste0(plot_dir, sample_id, "_", cluster_of_interest, "_prop.pdf")
  
  # Remove existing file if it exists
  if (file.exists(plot_filename)) {
    file.remove(plot_filename)
  }
  
  # Save the new plot
  ggsave(filename = plot_filename, plot = my_plot, width = 10, height = 8)
} else {
  warning(paste("No data found for", cluster_of_interest, "in", sample_id))
}

# Write results to CSV
csv_filename <- paste0(file_dir, sample_id, "_cluster_props.csv")
write.csv(results_df, csv_filename, row.names = FALSE)

# --------------------------------------------------------
# **Step 10: Correlation Analysis Between Clusters**
# --------------------------------------------------------

# Filter out rows with NA or empty clusters
results_df_filtered <- results_df %>%
  filter(!is.na(cluster) & cluster != "")

# Pivot the data to wide format
wide_results <- results_df_filtered %>%
  select(panel_id, cluster, proportion) %>%
  pivot_wider(names_from = cluster, values_from = proportion, values_fill = 0)

# Remove panels with zero assigned cells (if any)
wide_results <- wide_results %>%
  filter(rowSums(select(., -panel_id)) > 0)

# Check if there are at least two clusters
cluster_names <- colnames(wide_results)[-1]  # Exclude 'panel_id'

# Check if there are at least two clusters
print("Number of clusters found:")
print(length(cluster_names))
if (length(cluster_names) < 2) {
  stop("Not enough clusters to compute correlations")
}

# ----------------------------------------------
# **Step 11: Calculate Correlation Matrix**
# ----------------------------------------------

# Calculate correlation matrix
cor_matrix <- cor(wide_results[,-1])  # Exclude 'panel_id' column
print("Correlation matrix dimensions:")
print(dim(cor_matrix))

# ----------------------------------------------
# **Step 12: Permutation Testing**
# ----------------------------------------------

set.seed(123)  # For reproducibility
n_permutations <- 1000
perm_results <- list()

# Run permutation tests
for (i in 1:(length(cluster_names) - 1)) {
  for (j in (i + 1):length(cluster_names)) {
    cluster_a <- cluster_names[i]
    cluster_b <- cluster_names[j]
    
    observed_corr <- cor(wide_results[[cluster_a]], wide_results[[cluster_b]])
    permuted_corrs <- numeric(n_permutations)
    
    for (k in 1:n_permutations) {
      shuffled_b <- sample(wide_results[[cluster_b]])
      permuted_corrs[k] <- cor(wide_results[[cluster_a]], shuffled_b)
    }
    
    p_value <- if (observed_corr >= 0) {
      sum(permuted_corrs >= observed_corr) / n_permutations
    } else {
      sum(permuted_corrs <= observed_corr) / n_permutations
    }
    
    perm_results[[paste(cluster_a, cluster_b, sep = "_vs_")]] <- data.frame(
      cluster_a = cluster_a,
      cluster_b = cluster_b,
      observed_corr = observed_corr,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }
}

# ----------------------------------------------
# **Step 13: Process Results**
# ----------------------------------------------

# Combine results and adjust p-values
perm_results_df <- bind_rows(perm_results)
perm_results_df$adj_p_value <- p.adjust(perm_results_df$p_value, method = "BH")

# Filter significant pairs
significant_pairs <- perm_results_df %>%
  filter(adj_p_value < p_value_threshold)

print("Number of significant pairs found:")
print(nrow(significant_pairs))

# Save results
signif_pairs_filename <- paste0(file_dir, sample_id, "_signif_pairs_pval", p_value_threshold * 100, ".csv")
write.csv(significant_pairs, signif_pairs_filename, row.names = FALSE)

# ----------------------------------------------
# **Step 14: Correlation Heatmap**
# ----------------------------------------------

# Create heatmap
cor_melted <- melt(cor_matrix)
heatmap_plot <- ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limit = c(-1, 1), name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  labs(title = paste("Cluster Co-localization Heatmap in", sample_id,
                     "(Panel:", panel_size, "um,", overlap * 100, "% overlap)"),
       x = "Cluster", y = "Cluster") +
  coord_fixed()

# Save heatmap
heatmap_filename <- paste0(plot_dir, sample_id, "_corr_heatmap.pdf")

if (file.exists(heatmap_filename)) {
  file.remove(heatmap_filename)
}

# Save the heatmap
ggsave(filename = heatmap_filename, plot = heatmap_plot, width = 10, height = 8)
# ----------------------------------------------
# **Step 15: Network Graph**
# ----------------------------------------------

if (nrow(significant_pairs) > 0) {
  # Prepare edges data
  edges <- significant_pairs %>%
    select(cluster_a, cluster_b, observed_corr)
  
  edges_fixed <- data.frame(
    cluster_a = edges$cluster_a,
    cluster_b = edges$cluster_b,
    weight = abs(edges$observed_corr),
    sign = ifelse(edges$observed_corr > 0, "Positive", "Negative")
  )
  
  # Create and plot network
  g <- graph_from_data_frame(edges_fixed, directed = FALSE)
  E(g)$color <- ifelse(edges_fixed$sign == "Positive", "red", "blue")
  
  network_graph_filename <- paste0(plot_dir, sample_id, "_network_graph.pdf")
  pdf(network_graph_filename, width = 10, height = 8)
  plot(g,
       edge.width = E(g)$weight * 5,
       edge.color = E(g)$color,
       vertex.size = 30,
       vertex.label.cex = 1.2,
       layout = layout_with_fr(g),
       main = paste("Network of Significant Cluster Co-localizations in", sample_id,
                    "\n(Panel:", panel_size, "um,", overlap * 100, "% overlap, p <", p_value_threshold, ")"))
  dev.off()
}

# ----------------------------------------------
# **Step 16: Spatial Distribution Plots**
# ----------------------------------------------

if (nrow(significant_pairs) > 0) {
  for (k in 1:nrow(significant_pairs)) {
    cluster_pair <- significant_pairs[k, c("cluster_a", "cluster_b")]
    
    coords_pair <- coords %>%
      filter(cluster %in% c(cluster_pair$cluster_a, cluster_pair$cluster_b))
    
    spatial_plot <- ggplot(coords_pair, aes(x = x, y = y, color = cluster)) +
      geom_point(alpha = 0.6, size = 0.5) +
      labs(title = paste("Spatial Distribution of",
                         cluster_pair$cluster_a, "and",
                         cluster_pair$cluster_b, "in", sample_id,
                         "\n(Panel:", panel_size, "um,", overlap * 100, "% overlap, p <", p_value_threshold, ")"),
           x = "X Coordinate", y = "Y Coordinate") +
      theme_minimal()
    
    spatial_plot_filename <- paste0(spatial_plot_dir, sample_id, "_spatial_",
                                    cluster_pair$cluster_a, "_",
                                    cluster_pair$cluster_b, ".pdf")
    
    # Remove existing file if it exists
    if (file.exists(spatial_plot_filename)) {
      file.remove(spatial_plot_filename)
    }
    
    # Save the spatial plot
    ggsave(filename = spatial_plot_filename, 
           plot = spatial_plot, 
           width = 10, 
           height = 8)
  }
}
