# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)

# Step 1: Load the Seurat object
JY10 <- readRDS("rds_files/JY10_object.rds")

# Step 2: Load the relabeled barcodes
relabeled_barcodes <- read_csv("relabeled_barcodes.csv")

# Step 3: Prepare the barcodes for matching

# Get cell names from Seurat object
seurat_cell_names <- colnames(JY10)

# Check if the 'cell' column exists in relabeled_barcodes
if (!"cell" %in% colnames(relabeled_barcodes)) {
  stop("The 'cell' column is not present in 'relabeled_barcodes'. Please check the column names.")
}

# Check how many barcodes match directly
direct_matches <- sum(relabeled_barcodes$cell %in% seurat_cell_names)
cat("Number of direct barcode matches:", direct_matches, "\n")

if (direct_matches == 0) {
  # No direct matches found, likely due to prefixes/suffixes
  # Investigate the format of cell names in Seurat object
  head(seurat_cell_names)
  
  # Example: If Seurat cell names have a prefix like "slice1.008um_"
  # Extract the prefix from Seurat cell names
  prefix <- sub("(.+_).+", "\\1", seurat_cell_names[1])
  cat("Detected prefix in Seurat cell names:", prefix, "\n")
  
  # Add the prefix to the barcodes in relabeled_barcodes
  relabeled_barcodes <- relabeled_barcodes %>%
    mutate(cell = paste0(prefix, cell))
  
  # Check matches again
  updated_matches <- sum(relabeled_barcodes$cell %in% seurat_cell_names)
  cat("Number of barcode matches after adding prefix:", updated_matches, "\n")
  
  if (updated_matches == 0) {
    stop("Barcodes in the CSV file do not match cell names in the Seurat object even after adding prefix. Please check for discrepancies.")
  }
}

# Step 4: Add labels as metadata to the Seurat object

# Ensure that the cell barcodes are set as rownames
cell_labels <- relabeled_barcodes %>%
  dplyr::select(cell, Label) %>%
  column_to_rownames(var = "cell")

# Verify that all cells in cell_labels are in the Seurat object
missing_cells <- setdiff(rownames(cell_labels), seurat_cell_names)
if (length(missing_cells) > 0) {
  cat("The following cells are in 'relabeled_barcodes' but not in the Seurat object:\n")
  print(head(missing_cells))
  stop("Please ensure that the barcodes match between the CSV file and the Seurat object.")
}

# Add the labels to the Seurat object's metadata
JY10 <- AddMetaData(JY10, metadata = cell_labels)

# Step 5: Subset the Seurat object into tumor and non-tumor objects

# Verify that the 'Label' column has been added to the metadata
if (!"Label" %in% colnames(JY10@meta.data)) {
  stop("The 'Label' column was not successfully added to the Seurat object's metadata.")
}

# Create the tumor subset
JY10_tumor <- subset(JY10, subset = Label == "Tumor")

# Create the non-tumor subset
JY10_nontumor <- subset(JY10, subset = Label == "Non-Tumor")

# Step 6: Save the tumor and non-tumor Seurat objects
saveRDS(JY10_tumor, file = "JY10_tumor.rds")
saveRDS(JY10_nontumor, file = "JY10_nontumor.rds")

# Step 7: Output summary information
cat("Number of tumor cells:", ncol(JY10_tumor), "\n")
cat("Number of non-tumor cells:", ncol(JY10_nontumor), "\n")