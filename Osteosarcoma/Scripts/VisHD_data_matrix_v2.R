# Load necessary libraries
library(Seurat)      # For handling single-cell RNA-seq data
library(Matrix)      # For matrix operations
library(data.table)  # For efficient data manipulation

# Function to preprocess counts data (compatible with Seurat v5)
preprocess_counts_data <- function(seurat_object, assay_name = "Spatial.008um") {
  # Set default assay
  DefaultAssay(seurat_object) <- assay_name
  
  # Access the counts matrix from the counts slot
  counts_matrix <- GetAssayData(seurat_object, layer = "counts")
  
  # Transpose the counts matrix
  counts_matrix_t <- t(counts_matrix)
  
  # Convert the sparse matrix to a data frame of non-zero entries
  counts_df <- as.data.frame(summary(counts_matrix_t))
  
  # Add barcode and gene names based on row and column indices
  counts_df$Barcode <- rownames(counts_matrix_t)[counts_df$i]
  counts_df$Gene <- colnames(counts_matrix_t)[counts_df$j]
  
  # Select and rename relevant columns
  counts_df <- counts_df[, c("Barcode", "Gene", "x")]
  colnames(counts_df)[3] <- "Count"
  
  # Ensure consistent gene naming (e.g., uppercase)
  counts_df$Gene <- toupper(counts_df$Gene)
  
  return(counts_df)
}

# Function to extract cell types and their top marker genes from Markers file
get_cell_types <- function(Markers, p_val_threshold = 1e-5, top_n = NULL) {
  # Load necessary library
  library(data.table)
  
  # Read Markers data
  Markers_dt <- fread(Markers)
  
  # Ensure Markers has the expected columns
  required_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  if (!all(required_cols %in% colnames(Markers_dt))) {
    stop("Markers file must contain the following columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Ensure consistent gene naming (e.g., uppercase)
  Markers_dt$gene <- toupper(Markers_dt$gene)
  
  # Filter genes based on p-value threshold
  Markers_dt <- Markers_dt[p_val_adj < p_val_threshold]
  
  # Order Markers data by cluster and decreasing avg_log2FC
  Markers_dt <- Markers_dt[order(cluster, -avg_log2FC)]
  
  # For each cluster, get top_n genes by avg_log2FC
  if (!is.null(top_n)) {
    cell_types_dt <- Markers_dt[, head(.SD, top_n), by = cluster]
  } else {
    cell_types_dt <- Markers_dt
  }
  
  # Create a list where names are clusters and values are vectors of genes
  cell_types_list <- cell_types_dt[, .(genes = list(unique(gene))), by = cluster]
  
  # Convert to named list
  cell_types <- setNames(cell_types_list$genes, cell_types_list$cluster)
  
  return(cell_types)
}

# Function to normalize counts and compute cell type scores
normalize_and_score <- function(csv_file_path, Markers, output_file_path = NULL, p_val_threshold = 1e-5, top_n = NULL) {
  # Load necessary library
  library(data.table)
  
  # Step 1: Read the counts CSV file
  counts_dt <- fread(csv_file_path)
  
  # Ensure the data has the expected columns
  if (!all(c("Barcode", "Gene", "Count") %in% colnames(counts_dt))) {
    stop("Input CSV must contain 'Barcode', 'Gene', and 'Count' columns.")
  }
  
  # Ensure consistent gene naming (e.g., uppercase)
  counts_dt$Gene <- toupper(counts_dt$Gene)
  
  # Step 2: Compute min and max per gene
  gene_stats <- counts_dt[, .(Min = min(Count), Max = max(Count)), by = Gene]
  
  # Step 3: Merge min and max back to counts_dt
  counts_dt <- merge(counts_dt, gene_stats, by = "Gene")
  
  # Step 4: Perform min-max normalization per gene
  counts_dt[, NormCount := ifelse(Max - Min == 0, 0, (Count - Min) / (Max - Min))]
  
  # Step 5: Define cell types and their marker genes using Markers
  cell_types <- get_cell_types(Markers, p_val_threshold, top_n)
  
  # Step 6: Initialize scores data.table with unique barcodes
  barcodes <- unique(counts_dt$Barcode)
  scores_dt <- data.table(Barcode = barcodes)
  
  # Step 7: Compute scores for each cell type
  for (cell_type in names(cell_types)) {
    genes <- cell_types[[cell_type]]
    # Filter counts_dt for the marker genes of the current cell type
    marker_counts <- counts_dt[Gene %in% genes, .(Score = sum(NormCount)), by = Barcode]
    # Merge the scores back to the scores_dt table
    scores_dt <- merge(scores_dt, marker_counts, by = "Barcode", all.x = TRUE)
    # Rename the 'Score' column to the current cell type
    setnames(scores_dt, "Score", cell_type)
  }
  
  # Step 8: Replace NA values with 0 (barcodes without marker genes for certain cell types)
  for (cell_type in names(cell_types)) {
    scores_dt[is.na(get(cell_type)), (cell_type) := 0]
  }
  
  # Optional: Write the scores to a CSV file
  if (!is.null(output_file_path)) {
    fwrite(scores_dt, file = output_file_path)
  }
  
  return(scores_dt)
}

# Function to annotate barcodes based on highest cell type score
annotate_barcodes <- function(scores_dt, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Melt the scores data.table to long format
  scores_long <- melt(scores_dt, id.vars = "Barcode", variable.name = "CellType", value.name = "Score")
  
  # For each barcode, find the cell type(s) with the highest score
  annotations <- scores_long[, {
    max_score <- max(Score)
    if (max_score == 0) {
      # If all scores are zero, label as "unassigned"
      .(Graph_based = "unassigned")
    } else {
      top_cell_types <- CellType[Score == max_score]
      .(Graph_based = paste(top_cell_types, collapse = " / "))
    }
  }, by = Barcode]
  
  # Optional: Write to CSV
  if (!is.null(output_file_path)) {
    fwrite(annotations, file = output_file_path)
  }
  
  return(annotations)
}

# Read the Seurat objects
JY4_data <- readRDS("rds_files/JY4_object.rds")
JY8_data <- readRDS("rds_files/JY8_object.rds")
JY10_data <- readRDS("rds_files/JY10_object.rds")
JY23_data <- readRDS("rds_files/JY23_object.rds")

# Preprocess counts data
JY4_counts_df <- preprocess_counts_data(JY4_data, assay_name = "Spatial.008um")
JY8_counts_df <- preprocess_counts_data(JY8_data, assay_name = "Spatial.008um")
JY10_counts_df <- preprocess_counts_data(JY10_data, assay_name = "Spatial.008um")
JY23_counts_df <- preprocess_counts_data(JY23_data, assay_name = "Spatial.008um")

# Write the data frames to CSV files
fwrite(JY4_counts_df, file = "JY4_counts_matrix.csv")
fwrite(JY8_counts_df, file = "JY8_counts_matrix.csv")
fwrite(JY10_counts_df, file = "JY10_counts_matrix.csv")
fwrite(JY23_counts_df, file = "JY23_counts_matrix.csv")

# Set the path to your marker genes CSV file (e.g., Markers)
Markers <- "Deduplicated_Gene_Data.csv"  # Replace with the actual path to your marker genes CSV file

# Set the p-value threshold for selecting marker genes
p_val_threshold <- 1e-5  # Adjust as needed for stringency

# Optionally set top N genes per cell type
top_n_genes <- 100  # Adjust as needed; set to NULL to include all genes passing the p-value threshold

# Run the normalize_and_score function for each dataset
JY4_scores <- normalize_and_score("JY4_counts_matrix.csv", Markers, "JY4_cell_type_scores.csv", p_val_threshold, top_n_genes)
JY8_scores <- normalize_and_score("JY8_counts_matrix.csv", Markers, "JY8_cell_type_scores.csv", p_val_threshold, top_n_genes)
JY10_scores <- normalize_and_score("JY10_counts_matrix.csv", Markers, "JY10_cell_type_scores.csv", p_val_threshold, top_n_genes)
JY23_scores <- normalize_and_score("JY23_counts_matrix.csv", Markers, "JY23_cell_type_scores.csv", p_val_threshold, top_n_genes)

# Annotate barcodes for each dataset
JY4_annotations <- annotate_barcodes(JY4_scores, "JY4_annotations.csv")
JY8_annotations <- annotate_barcodes(JY8_scores, "JY8_annotations.csv")
JY10_annotations <- annotate_barcodes(JY10_scores, "JY10_annotations.csv")
JY23_annotations <- annotate_barcodes(JY23_scores, "JY23_annotations.csv")

# Function to reorder annotations based on barcodes order
reorder_annotations <- function(annotations_file_path, barcodes_file_path, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Read the barcodes file
  barcodes_dt <- fread(barcodes_file_path)
  
  # Read the annotations file
  annotations_dt <- fread(annotations_file_path)
  
  # Ensure that both data.tables have the 'Barcode' column
  if (!"Barcode" %in% colnames(barcodes_dt)) {
    stop("Barcodes file must contain a 'Barcode' column.")
  }
  
  if (!"Barcode" %in% colnames(annotations_dt)) {
    stop("Annotations file must contain a 'Barcode' column.")
  }
  
  # Merge annotations into barcodes, preserving the order of barcodes
  merged_dt <- merge(barcodes_dt, annotations_dt, by = "Barcode", all.x = TRUE, sort = FALSE)
  
  # Optional: Write to CSV
  if (!is.null(output_file_path)) {
    fwrite(merged_dt, file = output_file_path)
  }
  
  return(merged_dt)
}

# Reorder annotations for each dataset
JY4_reordered_annotations <- reorder_annotations("JY4_annotations.csv", "JY4_barcodes.csv", "JY4_annotations_reordered.csv")
JY8_reordered_annotations <- reorder_annotations("JY8_annotations.csv", "JY8_barcodes.csv", "JY8_annotations_reordered.csv")
JY10_reordered_annotations <- reorder_annotations("JY10_annotations.csv", "JY10_barcodes.csv", "JY10_annotations_reordered.csv")
JY23_reordered_annotations <- reorder_annotations("JY23_annotations.csv", "JY23_barcodes.csv", "JY23_annotations_reordered.csv")

# View the first few rows of the annotations
head(JY4_annotations)

# View the counts of barcodes annotated for each cell type
count_cell_types <- function(annotations_dt) {
  # Copy the annotations data.table to avoid modifying the original
  annotations_dt <- copy(annotations_dt)
  
  # Split 'Graph_based' into individual cell types for each barcode
  annotations_dt[, CellTypes := strsplit(Graph_based, " / ")]
  
  # Unnest the 'CellTypes' list column to have one cell type per row
  annotations_long <- annotations_dt[, .(Barcode, CellType = unlist(CellTypes))]
  
  # Count the number of barcodes per cell type
  cell_type_counts <- annotations_long[, .N, by = CellType]
  
  # Order the results by the number of barcodes in descending order
  cell_type_counts <- cell_type_counts[order(-N)]
  
  return(cell_type_counts)
}

# Apply the function to each dataset
JY4_cell_type_counts <- count_cell_types(JY4_annotations)
print(JY4_cell_type_counts)

JY8_cell_type_counts <- count_cell_types(JY8_annotations)
print(JY8_cell_type_counts)

JY10_cell_type_counts <- count_cell_types(JY10_annotations)
print(JY10_cell_type_counts)

JY23_cell_type_counts <- count_cell_types(JY23_annotations)
print(JY23_cell_type_counts)

################################## Updated scoring 
# Load necessary libraries
library(Seurat)      # For handling single-cell RNA-seq data
library(Matrix)      # For matrix operations
library(data.table)  # For efficient data manipulation

# Function to preprocess counts data (compatible with Seurat v5)
preprocess_counts_data <- function(seurat_object, assay_name = "Spatial.008um") {
  # Set default assay
  DefaultAssay(seurat_object) <- assay_name
  
  # Access the counts matrix from the counts slot
  counts_matrix <- GetAssayData(seurat_object, layer = "counts")
  
  # Transpose the counts matrix
  counts_matrix_t <- t(counts_matrix)
  
  # Convert the sparse matrix to a data frame of non-zero entries
  counts_df <- as.data.frame(summary(counts_matrix_t))
  
  # Add barcode and gene names based on row and column indices
  counts_df$Barcode <- rownames(counts_matrix_t)[counts_df$i]
  counts_df$Gene <- colnames(counts_matrix_t)[counts_df$j]
  
  # Select and rename relevant columns
  counts_df <- counts_df[, c("Barcode", "Gene", "x")]
  colnames(counts_df)[3] <- "Count"
  
  # Ensure consistent gene naming (e.g., uppercase)
  counts_df$Gene <- toupper(counts_df$Gene)
  
  return(counts_df)
}

# Function to extract cell types and their marker genes from marker file and compute combined weights
get_cell_types_with_combined_weights <- function(marker_file_path, p_val_threshold = 1e-5, top_n = NULL) {
  # Load necessary library
  library(data.table)
  
  # Read marker data
  marker_dt <- fread(marker_file_path)
  
  # Ensure the marker file has the required columns
  required_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
  if (!all(required_cols %in% colnames(marker_dt))) {
    stop("Marker file must contain the following columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Ensure consistent gene naming (e.g., uppercase)
  marker_dt$gene <- toupper(marker_dt$gene)
  
  # Filter genes based on p-value threshold
  marker_dt <- marker_dt[p_val_adj < p_val_threshold]
  
  # Order marker data by cluster and decreasing avg_log2FC
  marker_dt <- marker_dt[order(cluster, -avg_log2FC)]
  
  # For each cluster, get top_n genes by avg_log2FC
  if (!is.null(top_n)) {
    marker_dt <- marker_dt[, head(.SD, top_n), by = cluster]
  }
  
  # Calculate the number of cell types each gene is associated with (specificity weight)
  gene_counts <- marker_dt[, .N, by = gene]
  setnames(gene_counts, "N", "cell_type_count")
  
  # Merge the counts back to marker_dt
  marker_dt <- merge(marker_dt, gene_counts, by = "gene")
  
  # Calculate specificity weight as inverse of cell type count
  marker_dt[, specificity_weight := 1 / cell_type_count]
  
  # Normalize avg_log2FC across all genes
  max_logFC <- max(marker_dt$avg_log2FC, na.rm = TRUE)
  min_logFC <- min(marker_dt$avg_log2FC, na.rm = TRUE)
  marker_dt[, logFC_norm := ifelse(max_logFC - min_logFC == 0, 0, (avg_log2FC - min_logFC) / (max_logFC - min_logFC))]
  
  # Calculate combined weight
  marker_dt[, combined_weight := specificity_weight * logFC_norm]
  
  # Prepare data for splitting
  marker_data_to_split <- marker_dt[, .(cluster, gene, combined_weight)]
  
  # Split the data.table by 'cluster' and ensure each element is a data.table
  marker_list <- split(marker_data_to_split[, .(gene, combined_weight)], f = marker_data_to_split$cluster)
  marker_list <- lapply(marker_list, as.data.table)
  
  return(marker_list)
}


# Function to normalize counts and compute weighted cell type scores using combined weights
normalize_and_score_with_combined_weights <- function(csv_file_path, marker_list, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Step 1: Read the counts CSV file
  counts_dt <- fread(csv_file_path)
  
  # Ensure the data has the expected columns
  if (!all(c("Barcode", "Gene", "Count") %in% colnames(counts_dt))) {
    stop("Input CSV must contain 'Barcode', 'Gene', and 'Count' columns.")
  }
  
  # Ensure consistent gene naming (e.g., uppercase)
  counts_dt$Gene <- toupper(counts_dt$Gene)
  
  # Step 2: Compute min and max per gene
  gene_stats <- counts_dt[, .(Min = min(Count), Max = max(Count)), by = Gene]
  
  # Step 3: Merge min and max back to counts_dt
  counts_dt <- merge(counts_dt, gene_stats, by = "Gene")
  
  # Step 4: Perform min-max normalization per gene
  counts_dt[, NormCount := ifelse(Max - Min == 0, 0, (Count - Min) / (Max - Min))]
  
  # Step 5: Initialize scores data.table with unique barcodes
  barcodes <- unique(counts_dt$Barcode)
  scores_dt <- data.table(Barcode = barcodes)
  
  # Step 6: Compute weighted scores for each cell type using combined weights
  for (cell_type in names(marker_list)) {
    markers <- marker_list[[cell_type]]
    # Ensure markers is a data.table
    markers <- as.data.table(markers)
    # Merge counts with marker combined weights
    marker_counts <- counts_dt[Gene %in% markers$gene]
    if (nrow(marker_counts) == 0) {
      # No marker genes present for this cell type in counts data
      scores_dt[, (cell_type) := 0]
      next
    }
    marker_counts <- merge(marker_counts, markers, by.x = "Gene", by.y = "gene", all.x = FALSE)
    # Compute weighted NormCount using combined_weight
    marker_counts[, WeightedCount := NormCount * combined_weight]
    # Sum weighted counts per barcode
    marker_scores <- marker_counts[, .(Score = sum(WeightedCount)), by = Barcode]
    # Merge the scores back to the scores_dt table
    scores_dt <- merge(scores_dt, marker_scores, by = "Barcode", all.x = TRUE)
    # Rename the 'Score' column to the current cell type
    setnames(scores_dt, "Score", cell_type)
  }
  
  # Step 7: Replace NA values with 0 (barcodes without marker genes for certain cell types)
  for (cell_type in names(marker_list)) {
    scores_dt[is.na(get(cell_type)), (cell_type) := 0]
  }
  
  # Optional: Write the scores to a CSV file
  if (!is.null(output_file_path)) {
    fwrite(scores_dt, file = output_file_path)
  }
  
  return(scores_dt)
}

# Function to annotate barcodes based on highest cell type score
annotate_barcodes <- function(scores_dt, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Melt the scores data.table to long format
  scores_long <- melt(scores_dt, id.vars = "Barcode", variable.name = "CellType", value.name = "Score")
  
  # For each barcode, find the cell type(s) with the highest score
  annotations <- scores_long[, {
    max_score <- max(Score)
    if (max_score == 0) {
      # If all scores are zero, label as "unassigned"
      .(Graph_based = "unassigned")
    } else {
      top_cell_types <- CellType[Score == max_score]
      .(Graph_based = paste(top_cell_types, collapse = " / "))
    }
  }, by = Barcode]
  
  # Optional: Write to CSV
  if (!is.null(output_file_path)) {
    fwrite(annotations, file = output_file_path)
  }
  
  return(annotations)
}

# Function to reorder annotations based on barcodes order
reorder_annotations <- function(annotations_file_path, barcodes_file_path, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Read the barcodes file
  barcodes_dt <- fread(barcodes_file_path)
  
  # Read the annotations file
  annotations_dt <- fread(annotations_file_path)
  
  # Ensure that both data.tables have the 'Barcode' column
  if (!"Barcode" %in% colnames(barcodes_dt)) {
    stop("Barcodes file must contain a 'Barcode' column.")
  }
  
  if (!"Barcode" %in% colnames(annotations_dt)) {
    stop("Annotations file must contain a 'Barcode' column.")
  }
  
  # Merge annotations into barcodes, preserving the order of barcodes
  merged_dt <- merge(barcodes_dt, annotations_dt, by = "Barcode", all.x = TRUE, sort = FALSE)
  
  # Optional: Write to CSV
  if (!is.null(output_file_path)) {
    fwrite(merged_dt, file = output_file_path)
  }
  
  return(merged_dt)
}

# Function to count the number of barcodes annotated for each cell type
count_cell_types <- function(annotations_dt) {
  # Copy the annotations data.table to avoid modifying the original
  annotations_dt <- copy(annotations_dt)
  
  # Split 'Graph_based' into individual cell types for each barcode
  annotations_dt[, CellTypes := strsplit(Graph_based, " / ")]
  
  # Unnest the 'CellTypes' list column to have one cell type per row
  annotations_long <- annotations_dt[, .(Barcode, CellType = unlist(CellTypes))]
  
  # Count the number of barcodes per cell type
  cell_type_counts <- annotations_long[, .N, by = CellType]
  
  # Order the results by the number of barcodes in descending order
  cell_type_counts <- cell_type_counts[order(-N)]
  
  return(cell_type_counts)
}

# Read the Seurat objects
JY4_data <- readRDS("rds_files/JY4_object.rds")
JY8_data <- readRDS("rds_files/JY8_object.rds")
JY10_data <- readRDS("rds_files/JY10_object.rds")
JY23_data <- readRDS("rds_files/JY23_object.rds")

# Preprocess counts data
JY4_counts_df <- preprocess_counts_data(JY4_data, assay_name = "Spatial.008um")
JY8_counts_df <- preprocess_counts_data(JY8_data, assay_name = "Spatial.008um")
JY10_counts_df <- preprocess_counts_data(JY10_data, assay_name = "Spatial.008um")
JY23_counts_df <- preprocess_counts_data(JY23_data, assay_name = "Spatial.008um")

# Write the data frames to CSV files
fwrite(JY4_counts_df, file = "JY4_counts_matrix.csv")
fwrite(JY8_counts_df, file = "JY8_counts_matrix.csv")
fwrite(JY10_counts_df, file = "JY10_counts_matrix.csv")
fwrite(JY23_counts_df, file = "JY23_counts_matrix.csv")

# Set the path to your marker genes CSV file
Markers <- "Deduplicated_Gene_Data.csv"  # Replace with the actual path

# Set the p-value threshold for selecting marker genes
p_val_threshold <- 1e-5  # Adjust as needed

# Optionally set top N genes per cell type
top_n_genes <- NULL  # Adjust as needed; set to NULL to include all genes passing the p-value threshold

# Load the marker list with weights using the corrected function
marker_list <- get_cell_types_with_weights(Markers, p_val_threshold, top_n_genes)

# Run the normalize_and_score_with_weights function for each dataset
JY4_scores <- normalize_and_score_with_weights("JY4_counts_matrix.csv", marker_list, "JY4_cell_type_scores.csv")
JY8_scores <- normalize_and_score_with_weights("JY8_counts_matrix.csv", marker_list, "JY8_cell_type_scores.csv")
JY10_scores <- normalize_and_score_with_weights("JY10_counts_matrix.csv", marker_list, "JY10_cell_type_scores.csv")
JY23_scores <- normalize_and_score_with_weights("JY23_counts_matrix.csv", marker_list, "JY23_cell_type_scores.csv")

# Annotate barcodes for each dataset
JY4_annotations <- annotate_barcodes(JY4_scores, "JY4_annotations.csv")
JY8_annotations <- annotate_barcodes(JY8_scores, "JY8_annotations.csv")
JY10_annotations <- annotate_barcodes(JY10_scores, "JY10_annotations.csv")
JY23_annotations <- annotate_barcodes(JY23_scores, "JY23_annotations.csv")

# Reorder annotations for each dataset
JY4_reordered_annotations <- reorder_annotations("JY4_annotations.csv", "JY4_barcodes.csv", "JY4_annotations_reordered.csv")
JY8_reordered_annotations <- reorder_annotations("JY8_annotations.csv", "JY8_barcodes.csv", "JY8_annotations_reordered.csv")
JY10_reordered_annotations <- reorder_annotations("JY10_annotations.csv", "JY10_barcodes.csv", "JY10_annotations_reordered.csv")
JY23_reordered_annotations <- reorder_annotations("JY23_annotations.csv", "JY23_barcodes.csv", "JY23_annotations_reordered.csv")

# View the first few rows of the annotations
head(JY4_annotations)

# Apply the function to each dataset
JY4_cell_type_counts <- count_cell_types(JY4_annotations)
print(JY4_cell_type_counts)

JY8_cell_type_counts <- count_cell_types(JY8_annotations)
print(JY8_cell_type_counts)

JY10_cell_type_counts <- count_cell_types(JY10_annotations)
print(JY10_cell_type_counts)

JY23_cell_type_counts <- count_cell_types(JY23_annotations)
print(JY23_cell_type_counts)

# Function to add annotations to Seurat object
add_annotations <- function(seurat_obj, annotations_dt) {
  # Ensure that barcodes in annotations match those in the Seurat object
  barcodes_in_obj <- colnames(seurat_obj)
  
  # Adjust barcodes in annotations to match Seurat object barcodes
  annotations_dt$Barcode <- as.character(annotations_dt$Barcode)
  
  # Find barcodes common to both annotations and Seurat object
  common_barcodes <- intersect(annotations_dt$Barcode, barcodes_in_obj)
  
  # Filter annotations to include only common barcodes
  annotations_dt <- annotations_dt[Barcode %in% common_barcodes]
  
  # Set identities in the Seurat object
  Idents(seurat_obj, cells = annotations_dt$Barcode) <- annotations_dt$Graph_based
  
  return(seurat_obj)
}

# Add annotations to each Seurat object
JY4_data <- add_annotations(JY4_data, JY4_annotations)
JY8_data <- add_annotations(JY8_data, JY8_annotations)
JY10_data <- add_annotations(JY10_data, JY10_annotations)
JY23_data <- add_annotations(JY23_data, JY23_annotations)

# Function to rename multi-cell-type annotations as "Multi"
simplify_annotations <- function(annotations_dt, output_file_path = NULL) {
  # Load necessary library
  library(data.table)
  
  # Copy the annotations data.table to avoid modifying the original
  annotations_dt <- copy(annotations_dt)
  
  # Identify rows where 'Graph_based' contains multiple cell types (indicated by " / ")
  multi_cell_type_rows <- grepl(" / ", annotations_dt$Graph_based)
  
  # Replace 'Graph_based' entries with "Multi" where multiple cell types are present
  annotations_dt[multi_cell_type_rows, Graph_based := "Multi"]
  
  # Optional: Write the simplified annotations to a CSV file
  if (!is.null(output_file_path)) {
    fwrite(annotations_dt, file = output_file_path)
  }
  
  return(annotations_dt)
}

# Simplify annotations for each dataset
JY4_annotations_simplified <- simplify_annotations(JY4_annotations, "JY4_annotations_simplified.csv")
JY8_annotations_simplified <- simplify_annotations(JY8_annotations, "JY8_annotations_simplified.csv")
JY10_annotations_simplified <- simplify_annotations(JY10_annotations, "JY10_annotations_simplified.csv")
JY23_annotations_simplified <- simplify_annotations(JY23_annotations, "JY23_annotations_simplified.csv")

# Update Seurat objects with simplified annotations
JY4_data <- add_annotations(JY4_data, JY4_annotations_simplified)
JY8_data <- add_annotations(JY8_data, JY8_annotations_simplified)
JY10_data <- add_annotations(JY10_data, JY10_annotations_simplified)
JY23_data <- add_annotations(JY23_data, JY23_annotations_simplified)

# Count cell types in simplified annotations
JY4_cell_type_counts_simplified <- count_cell_types(JY4_annotations_simplified)
print(JY4_cell_type_counts_simplified)

JY8_cell_type_counts_simplified <- count_cell_types(JY8_annotations_simplified)
print(JY8_cell_type_counts_simplified)

JY10_cell_type_counts_simplified <- count_cell_types(JY10_annotations_simplified)
print(JY10_cell_type_counts_simplified)

JY23_cell_type_counts_simplified <- count_cell_types(JY23_annotations_simplified)
print(JY23_cell_type_counts_simplified)



