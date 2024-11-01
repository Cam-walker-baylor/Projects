# Load necessary libraries
library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(writexl)
library(patchwork)

# Set base directory
base_dir <- "/home/walker/data/projects/osteosarcoma/data/mouse"

# Function to create and preprocess Seurat object for a single sample
create_and_preprocess_seurat_object <- function(data_dir, sample_name) {
  counts <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = counts, assay = "RNA")
  seurat_obj$sample <- sample_name
  
  # Use SCTransform for normalization
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  
  return(seurat_obj)
}

# Function to process a single mouse
process_mouse <- function(mouse_id) {
  mouse_dir <- file.path(base_dir, mouse_id)
  lung_dir <- file.path(mouse_dir, paste0(mouse_id, "L"))
  primary_dir <- file.path(mouse_dir, paste0(mouse_id, "P"))
  
  # Create and preprocess Seurat objects for lung and primary samples
  lung_seurat <- create_and_preprocess_seurat_object(lung_dir, "lung")
  primary_seurat <- create_and_preprocess_seurat_object(primary_dir, "primary")
  
  # List of Seurat objects
  seurat_list <- list(lung = lung_seurat, primary = primary_seurat)
  
  # Select integration features
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
  
  # Prepare the SCT assays for integration
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
                                    anchor.features = features)
  
  # Integrate data
  integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  
  # Proceed with downstream analysis on the integrated data
  DefaultAssay(integrated_data) <- "integrated"
  
  # Run PCA
  integrated_data <- RunPCA(integrated_data, verbose = FALSE)
  
  # Run UMAP and clustering
  integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)
  integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
  integrated_data <- FindClusters(integrated_data, resolution = 0.8)
  
  return(integrated_data)
}

# Process each mouse and collect integrated Seurat objects
mouse_ids <- c("M1", "M2", "M3")
seurat_objects <- list()
for (mouse in mouse_ids) {
  print(paste("Processing mouse", mouse))
  seurat_objects[[mouse]] <- process_mouse(mouse)
}

# Save integrated Seurat objects
output_dir <- "/data/u/walker/projects/osteosarcoma/data/mouse/seuratV4_rds_files" 
for (mouse in mouse_ids) {
  file_name <- file.path(output_dir, paste0(mouse, "_integrated_seurat_object.rds"))
  saveRDS(seurat_objects[[mouse]], file = file_name)
  print(paste("Saved", mouse, "Seurat object to", file_name))
}

# Load integrated Seurat objects
M1 <- readRDS(file.path(output_dir, "M1_integrated_seurat_object.rds"))
M2 <- readRDS(file.path(output_dir, "M2_integrated_seurat_object.rds"))
M3 <- readRDS(file.path(output_dir, "M3_integrated_seurat_object.rds"))

# Read the marker gene catalog CSV file
marker_gene_catalog <- read.csv("JY4_marker_catalog_with_Log2FC.csv")

# Function to create gene sets for a Seurat object based on the top 10 marker genes present in the data
create_gene_sets <- function(seurat_obj, marker_gene_catalog) {
  genes_in_data <- rownames(seurat_obj)
  gene_sets <- list()
  
  # Check if marker_gene_catalog has a column for marker importance
  if (!"avg_log2FC" %in% colnames(marker_gene_catalog)) {
    stop("The marker gene catalog must contain a column named 'avg_log2FC' indicating marker importance.")
  }
  
  # Get unique cell types
  cell_types <- unique(marker_gene_catalog$cluster)
  
  # Initialize a list to store the number of genes available for each cell type
  num_genes_available <- c()
  
  # First, determine the number of genes available for each cell type
  for (cell_type in cell_types) {
    # Get the marker genes for the current cell type
    marker_genes <- marker_gene_catalog %>%
      filter(cluster == cell_type) %>%
      arrange(desc(avg_log2FC))
    
    # Keep only genes present in the data
    genes_present <- intersect(marker_genes$gene, genes_in_data)
    
    # Store the number of genes available
    num_genes_available[cell_type] <- length(genes_present)
  }
  
  # Determine the minimum number of genes available across all cell types
  min_genes <- min(num_genes_available)
  
  # Ensure that min_genes is at least 2, since AddModuleScore requires at least two genes
  if (min_genes < 2) {
    stop("Not enough marker genes found in the data to calculate module scores. Minimum required is 2.")
  }
  
  # Now, create gene sets using the top 'min_genes' genes for each cell type
  for (cell_type in cell_types) {
    # Get the marker genes for the current cell type
    marker_genes <- marker_gene_catalog %>%
      filter(cluster == cell_type) %>%
      arrange(desc(avg_log2FC))
    
    # Keep only genes present in the data
    marker_genes_present <- marker_genes %>%
      filter(gene %in% genes_in_data)
    
    # Take the top 'min_genes' genes
    top_genes <- head(marker_genes_present$gene, min_genes)
    
    gene_sets[[cell_type]] <- top_genes
  }
  
  return(gene_sets)
}


# Updated add_module_scores function
add_module_scores <- function(seurat_obj, gene_sets) {
  DefaultAssay(seurat_obj) <- "integrated"  # or "RNA" if you prefer
  
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    
    if (length(genes) >= 2) {  # AddModuleScore requires at least two genes
      seurat_obj <- AddModuleScore(
        object = seurat_obj,
        features = list(genes),
        name = set_name
      )
      cat(paste("Added module score for", set_name, "with", length(genes), "genes\n"))
    } else {
      cat(paste("Warning: Not enough genes from", set_name, "found in the dataset to calculate module score\n"))
    }
  }
  
  return(seurat_obj)
}

# Function to annotate clusters based on dominant module scores
annotate_clusters <- function(seurat_obj) {
  # Extract module scores and cluster assignments
  module_score_cols <- grep("1$", colnames(seurat_obj@meta.data), value = TRUE)
  
  # Remove the "1" suffix from column names
  module_scores <- seurat_obj@meta.data[, module_score_cols]
  colnames(module_scores) <- gsub("1$", "", module_score_cols)
  
  # Add cluster information
  module_scores$seurat_clusters <- seurat_obj$seurat_clusters
  
  # Determine the dominant cell type for each cell
  module_scores$dominant_cell_type <- colnames(module_scores[, -ncol(module_scores)])[max.col(module_scores[, -ncol(module_scores)], ties.method = "first")]
  
  # Find the most common cell type in each cluster
  cluster_annotations <- module_scores %>%
    group_by(seurat_clusters) %>%
    count(dominant_cell_type) %>%
    top_n(1, n) %>%
    ungroup() %>%
    select(seurat_clusters, dominant_cell_type)
  
  # Add annotations to the Seurat object
  seurat_obj$cluster_annotation <- cluster_annotations$dominant_cell_type[match(seurat_obj$seurat_clusters, cluster_annotations$seurat_clusters)]
  
  return(seurat_obj)
}

# Function to plot cluster annotations
plot_cluster_annotations <- function(seurat_obj, obj_name) {
  p <- DimPlot(seurat_obj, 
               reduction = "umap", 
               group.by = "cluster_annotation", 
               label = TRUE, 
               label.size = 3, 
               repel = TRUE) +
    ggtitle(paste(obj_name, "Cluster Annotations")) +
    theme(legend.position = "right")
  
  return(p)
}

# Process each Seurat object with the new gene sets and updated functions
seurat_objects <- list(M1 = M1, M2 = M2, M3 = M3)
for (obj_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[obj_name]]
  
  # Create gene sets for this Seurat object
  gene_sets <- create_gene_sets(seurat_obj, marker_gene_catalog)
  
  # Add module scores using the gene sets
  seurat_obj <- add_module_scores(seurat_obj, gene_sets)
  
  # Annotate clusters based on module scores
  seurat_obj <- annotate_clusters(seurat_obj)
  
  # Update the Seurat object in the list
  seurat_objects[[obj_name]] <- seurat_obj
  
  # Save the updated Seurat object
  saveRDS(seurat_obj, file = paste0(obj_name, "_with_module_scores.rds"))
  
  # Create and save visualization
  plot <- plot_cluster_annotations(seurat_obj, obj_name)
  ggsave(filename = paste0(obj_name, "_cluster_annotations.pdf"), plot = plot, width = 12, height = 10)
}

# Update M1, M2, M3 with the new Seurat objects
M1 <- seurat_objects$M1
M2 <- seurat_objects$M2
M3 <- seurat_objects$M3

# Optional: Proceed with further analysis, such as finding markers and visualizations

# Function to combine clusters with the same cell type
combine_clusters <- function(seurat_obj) {
  seurat_obj$combined_annotation <- factor(seurat_obj$cluster_annotation)
  seurat_obj <- SetIdent(seurat_obj, value = "combined_annotation")
  return(seurat_obj)
}

# Apply the function to each Seurat object
M1 <- combine_clusters(M1)
M2 <- combine_clusters(M2)
M3 <- combine_clusters(M3)

# Function to find markers and create visualizations
analyze_and_visualize <- function(seurat_obj, obj_name) {
  # Find markers
  all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Save markers to a CSV file
  write.csv(all_markers, file = paste0(obj_name, "_all_markers.csv"), row.names = FALSE)
  
  # Create UMAP plot
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "combined_annotation", 
                       label = TRUE, repel = TRUE) + 
    ggtitle(paste(obj_name, "Combined Cell Types")) +
    theme(legend.position = "right")
  
  # Create heatmap of top markers
  top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  heatmap_plot <- DoHeatmap(seurat_obj, features = top10$gene) + 
    ggtitle(paste(obj_name, "Top 10 Markers per Cell Type"))
  
  # Combine plots
  combined_plot <- umap_plot + heatmap_plot + plot_layout(ncol = 1)
  
  # Save the combined plot
  ggsave(paste0(obj_name, "_combined_plot.pdf"), combined_plot, width = 16, height = 20)
  
  # Return the markers and plots
  return(list(markers = all_markers, umap = umap_plot, heatmap = heatmap_plot))
}

# Analyze and visualize each Seurat object
M1_results <- analyze_and_visualize(M1, "M1")
M2_results <- analyze_and_visualize(M2, "M2")
M3_results <- analyze_and_visualize(M3, "M3")

# Display UMAP plots
print(M1_results$umap)
print(M2_results$umap)
print(M3_results$umap)

# Save the annotated Seurat objects
saveRDS(M1, file = "M1_annotated.rds")
saveRDS(M2, file = "M2_annotated.rds")
saveRDS(M3, file = "M3_annotated.rds")

cat("All processing complete!\n")

# Optional: Further analysis and visualization can be added as needed



##################################### ALT (Relative module score annotation)

# Function to find DEGs for all clusters
find_cluster_markers <- function(seurat_obj) {
  cluster_markers <- FindAllMarkers(
    object = seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  return(cluster_markers)
}

# Load necessary library for enrichment analysis
library(fisher.test)  # For Fisher's exact test

# Function to perform enrichment analysis
perform_enrichment_analysis <- function(cluster_markers, gene_sets, seurat_obj) {
  # Get background genes (all genes detected in the dataset)
  background_genes <- rownames(seurat_obj)
  
  # Initialize a data frame to store enrichment results
  enrichment_results <- data.frame(
    cluster = character(),
    cell_type = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get unique clusters
  clusters <- unique(cluster_markers$cluster)
  
  for (cluster_id in clusters) {
    # Get DEGs for this cluster
    degs <- cluster_markers %>%
      filter(cluster == cluster_id) %>%
      pull(gene)
    
    # For each cell type, perform Fisher's exact test
    for (cell_type in names(gene_sets)) {
      # Marker genes for this cell type
      marker_genes <- gene_sets[[cell_type]]
      
      # Create contingency table
      # Genes in cluster DEGs and in marker genes
      overlap_genes <- intersect(degs, marker_genes)
      num_overlap <- length(overlap_genes)
      
      # Genes in DEGs but not in marker genes
      num_degs_not_markers <- length(setdiff(degs, marker_genes))
      
      # Genes not in DEGs but in marker genes
      num_markers_not_degs <- length(setdiff(marker_genes, degs))
      
      # Genes not in DEGs and not in marker genes
      num_rest <- length(setdiff(background_genes, union(degs, marker_genes)))
      
      # Construct contingency table
      contingency_table <- matrix(
        c(
          num_overlap,
          num_markers_not_degs,
          num_degs_not_markers,
          num_rest
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      # Perform Fisher's exact test
      test_result <- fisher.test(contingency_table, alternative = "greater")
      
      # Store the result
      enrichment_results <- rbind(
        enrichment_results,
        data.frame(
          cluster = cluster_id,
          cell_type = cell_type,
          p_value = test_result$p.value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(enrichment_results)
}

# Function to annotate clusters based on enrichment results
annotate_clusters_based_on_enrichment <- function(enrichment_results) {
  # For each cluster, find the cell type with the lowest p-value
  cluster_annotations <- enrichment_results %>%
    group_by(cluster) %>%
    filter(p_value == min(p_value)) %>%
    ungroup()
  
  # Handle ties by selecting the cell type with the highest odds ratio or any other criteria
  # For simplicity, we'll select the first one in case of ties
  cluster_annotations <- cluster_annotations %>%
    group_by(cluster) %>%
    slice(1) %>%
    ungroup()
  
  return(cluster_annotations)
}

plot_cluster_annotations <- function(seurat_obj, obj_name) {
  p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "cluster_annotation",
    label = TRUE,
    label.size = 3,
    repel = TRUE
  ) +
    ggtitle(paste(obj_name, "Cluster Annotations")) +
    theme(legend.position = "right")
  
  return(p)
}

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Proceed with your pipeline up to the point of module scoring

# Process each Seurat object with the new gene sets and updated functions
seurat_objects <- list(M1 = M1, M2 = M2, M3 = M3)
for (obj_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[obj_name]]
  
  cat("\nProcessing", obj_name, "\n")
  
  # Create gene sets for this Seurat object
  gene_sets <- create_gene_sets(seurat_obj, marker_gene_catalog)
  
  # Add module scores using the gene sets
  seurat_obj <- add_module_scores(seurat_obj, gene_sets)
  
  # Update the Seurat object in the list
  seurat_objects[[obj_name]] <- seurat_obj
  
  # Find cluster markers
  cluster_markers <- find_cluster_markers(seurat_obj)
  
  # Perform enrichment analysis
  enrichment_results <- perform_enrichment_analysis(cluster_markers, gene_sets, seurat_obj)
  
  # Annotate clusters based on enrichment
  cluster_annotations <- annotate_clusters_based_on_enrichment(enrichment_results)
  
  # Add annotations to the Seurat object
  seurat_obj$cluster_annotation <- cluster_annotations$cell_type[match(seurat_obj$seurat_clusters, cluster_annotations$cluster)]
  
  # Update the Seurat object in the list
  seurat_objects[[obj_name]] <- seurat_obj
  
  # Save the updated Seurat object
  saveRDS(seurat_obj, file = paste0(obj_name, "_with_enrichment_annotations.rds"))
  
  # Create and save visualization
  plot <- plot_cluster_annotations(seurat_obj, obj_name)
  ggsave(filename = paste0(obj_name, "_cluster_annotations.pdf"), plot = plot, width = 12, height = 10)
}
