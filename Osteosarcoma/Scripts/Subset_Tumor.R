# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)
library(tibble)


# Step 1: Load the Seurat object
JY4 <- readRDS("rds_files/JY4_object.rds")
JY8 <- readRDS("rds_files/JY8_object.rds")
JY10 <- readRDS("rds_files/JY10_object.rds")
JY23 <- readRDS("rds_files/JY23_object.rds")


# Step 2: Load the relabeled barcodes
JY4_relabeled_barcodes <- read_csv("JY4_relabeled_barcodes.csv")
JY8_relabeled_barcodes <- read_csv("JY8_relabeled_barcodes.csv")
JY10_relabeled_barcodes <- read_csv("JY10_relabeled_barcodes.csv")
JY23_relabeled_barcodes <- read_csv("JY23_relabeled_barcodes.csv")


# Step 3: Prepare the barcodes for matching

# Get cell names from Seurat object
JY4_seurat_cell_names <- colnames(JY4)
JY8_seurat_cell_names <- colnames(JY8)
JY10_seurat_cell_names <- colnames(JY10)
JY23_seurat_cell_names <- colnames(JY23)


# Check if the 'cell' column exists in relabeled_barcodes
if (!"cell" %in% colnames(JY4_relabeled_barcodes)) {
  stop("The 'cell' column is not present in 'JY4_relabeled_barcodes'. Please check the column names.")
}

if (!"cell" %in% colnames(JY8_relabeled_barcodes)) {
  stop("The 'cell' column is not present in 'JY8_relabeled_barcodes'. Please check the column names.")
}

if (!"cell" %in% colnames(JY10_relabeled_barcodes)) {
  stop("The 'cell' column is not present in 'JY10_relabeled_barcodes'. Please check the column names.")
}

if (!"cell" %in% colnames(JY23_relabeled_barcodes)) {
  stop("The 'cell' column is not present in 'JY23_relabeled_barcodes'. Please check the column names.")
}

# Check how many barcodes match directly
JY4_direct_matches <- sum(JY4_relabeled_barcodes$cell %in% JY4_seurat_cell_names)
cat("Number of direct barcode matches:", JY4_direct_matches, "\n")

if (JY4_direct_matches == 0) {
  # No direct matches found, likely due to prefixes/suffixes
  # Investigate the format of cell names in Seurat object
  head(JY4_seurat_cell_names)
  
  # Example: If Seurat cell names have a prefix like "slice1.008um_"
  # Extract the prefix from Seurat cell names
  prefix <- sub("(.+_).+", "\\1", JY4_seurat_cell_names[1])
  cat("Detected prefix in Seurat cell names:", prefix, "\n")
  
  # Add the prefix to the barcodes in relabeled_barcodes
  JY4_relabeled_barcodes <- JY4_relabeled_barcodes %>%
    mutate(cell = paste0(prefix, cell))
  
  # Check matches again
  JY4_updated_matches <- sum(JY4_relabeled_barcodes$cell %in% JY4_seurat_cell_names)
  cat("Number of barcode matches after adding prefix:", JY4_updated_matches, "\n")
  
  if (JY4_updated_matches == 0) {
    stop("Barcodes in the CSV file do not match cell names in the Seurat object even after adding prefix. Please check for discrepancies.")
  }
}

JY8_direct_matches <- sum(JY8_relabeled_barcodes$cell %in% JY8_seurat_cell_names)
cat("Number of direct barcode matches:", JY8_direct_matches, "\n")

if (JY8_direct_matches == 0) {
  # No direct matches found, likely due to prefixes/suffixes
  # Investigate the format of cell names in Seurat object
  head(JY8_seurat_cell_names)
  
  # Example: If Seurat cell names have a prefix like "slice1.008um_"
  # Extract the prefix from Seurat cell names
  prefix <- sub("(.+_).+", "\\1", JY8_seurat_cell_names[1])
  cat("Detected prefix in Seurat cell names:", prefix, "\n")
  
  # Add the prefix to the barcodes in relabeled_barcodes
  JY8_relabeled_barcodes <- JY8_relabeled_barcodes %>%
    mutate(cell = paste0(prefix, cell))
  
  # Check matches again
  JY8_updated_matches <- sum(JY8_relabeled_barcodes$cell %in% JY8_seurat_cell_names)
  cat("Number of barcode matches after adding prefix:", JY8_updated_matches, "\n")
  
  if (JY8_updated_matches == 0) {
    stop("Barcodes in the CSV file do not match cell names in the Seurat object even after adding prefix. Please check for discrepancies.")
  }
}

JY10_direct_matches <- sum(JY10_relabeled_barcodes$cell %in% JY10_seurat_cell_names)
cat("Number of direct barcode matches:", JY10_direct_matches, "\n")

if (JY10_direct_matches == 0) {
  # No direct matches found, likely due to prefixes/suffixes
  # Investigate the format of cell names in Seurat object
  head(JY10_seurat_cell_names)
  
  # Example: If Seurat cell names have a prefix like "slice1.008um_"
  # Extract the prefix from Seurat cell names
  prefix <- sub("(.+_).+", "\\1", JY10_seurat_cell_names[1])
  cat("Detected prefix in Seurat cell names:", prefix, "\n")
  
  # Add the prefix to the barcodes in relabeled_barcodes
  JY10_relabeled_barcodes <- JY10_relabeled_barcodes %>%
    mutate(cell = paste0(prefix, cell))
  
  # Check matches again
  JY10_updated_matches <- sum(JY10_relabeled_barcodes$cell %in% JY10_seurat_cell_names)
  cat("Number of barcode matches after adding prefix:", JY10_updated_matches, "\n")
  
  if (JY10_updated_matches == 0) {
    stop("Barcodes in the CSV file do not match cell names in the Seurat object even after adding prefix. Please check for discrepancies.")
  }
}

JY23_direct_matches <- sum(JY23_relabeled_barcodes$cell %in% JY23_seurat_cell_names)
cat("Number of direct barcode matches:", JY23_direct_matches, "\n")

if (JY23_direct_matches == 0) {
  # No direct matches found, likely due to prefixes/suffixes
  # Investigate the format of cell names in Seurat object
  head(JY23_seurat_cell_names)
  
  # Example: If Seurat cell names have a prefix like "slice1.008um_"
  # Extract the prefix from Seurat cell names
  prefix <- sub("(.+_).+", "\\1", JY23_seurat_cell_names[1])
  cat("Detected prefix in Seurat cell names:", prefix, "\n")
  
  # Add the prefix to the barcodes in relabeled_barcodes
  JY23_relabeled_barcodes <- JY23_relabeled_barcodes %>%
    mutate(cell = paste0(prefix, cell))
  
  # Check matches again
  JY23_updated_matches <- sum(JY23_relabeled_barcodes$cell %in% JY23_seurat_cell_names)
  cat("Number of barcode matches after adding prefix:", JY23_updated_matches, "\n")
  
  if (JY23_updated_matches == 0) {
    stop("Barcodes in the CSV file do not match cell names in the Seurat object even after adding prefix. Please check for discrepancies.")
  }
}

# Step 4: Add labels as metadata to the Seurat object

# Ensure that the cell barcodes are set as rownames
JY4_cell_labels <- JY4_relabeled_barcodes %>%
  dplyr::select(cell, Label) %>%
  column_to_rownames(var = "cell")

JY8_cell_labels <- JY8_relabeled_barcodes %>%
  dplyr::select(cell, Label) %>%
  column_to_rownames(var = "cell")

JY10_cell_labels <- JY10_relabeled_barcodes %>%
  dplyr::select(cell, Label) %>%
  column_to_rownames(var = "cell")

JY23_cell_labels <- JY23_relabeled_barcodes %>%
  dplyr::select(cell, Label) %>%
  column_to_rownames(var = "cell")

# Verify that all cells in cell_labels are in the Seurat object
# JY4
missing_cells_JY4 <- setdiff(rownames(JY4_cell_labels), JY4_seurat_cell_names)
if (length(missing_cells_JY4) > 0) {
  cat("The following cells are in 'JY4_relabeled_barcodes' but not in the Seurat object:\n")
  print(head(missing_cells_JY4))
  stop("Please ensure that the barcodes match between the CSV file and the Seurat object.")
}

JY4 <- AddMetaData(JY4, metadata = JY4_cell_labels)

# JY8
missing_cells_JY8 <- setdiff(rownames(JY8_cell_labels), JY8_seurat_cell_names)
if (length(missing_cells_JY8) > 0) {
  cat("The following cells are in 'JY8_relabeled_barcodes' but not in the Seurat object:\n")
  print(head(missing_cells_JY8))
  stop("Please ensure that the barcodes match between the CSV file and the Seurat object.")
}

JY8 <- AddMetaData(JY8, metadata = JY8_cell_labels)

# JY10
missing_cells_JY10 <- setdiff(rownames(JY10_cell_labels), JY10_seurat_cell_names)
if (length(missing_cells_JY10) > 0) {
  cat("The following cells are in 'JY10_relabeled_barcodes' but not in the Seurat object:\n")
  print(head(missing_cells_JY10))
  stop("Please ensure that the barcodes match between the CSV file and the Seurat object.")
}

JY10 <- AddMetaData(JY10, metadata = JY10_cell_labels)

# JY23
missing_cells_JY23 <- setdiff(rownames(JY23_cell_labels), JY23_seurat_cell_names)
if (length(missing_cells_JY23) > 0) {
  cat("The following cells are in 'JY23_relabeled_barcodes' but not in the Seurat object:\n")
  print(head(missing_cells_JY23))
  stop("Please ensure that the barcodes match between the CSV file and the Seurat object.")
}

JY23 <- AddMetaData(JY23, metadata = JY23_cell_labels)

# Step 5: Subset the Seurat objects into tumor and non-tumor objects for all samples

# JY4
if (!"Label" %in% colnames(JY4@meta.data)) {
  stop("The 'Label' column was not successfully added to the JY4 Seurat object's metadata.")
}

JY4_tumor <- subset(JY4, subset = Label == "Tumor")
JY4_nontumor <- subset(JY4, subset = Label == "Non-Tumor")

# JY8
if (!"Label" %in% colnames(JY8@meta.data)) {
  stop("The 'Label' column was not successfully added to the JY8 Seurat object's metadata.")
}

JY8_tumor <- subset(JY8, subset = Label == "Tumor")
JY8_nontumor <- subset(JY8, subset = Label == "Non-Tumor")

# JY10
if (!"Label" %in% colnames(JY10@meta.data)) {
  stop("The 'Label' column was not successfully added to the JY10 Seurat object's metadata.")
}

JY10_tumor <- subset(JY10, subset = Label == "Tumor")
JY10_nontumor <- subset(JY10, subset = Label == "Non-Tumor")

# JY23
if (!"Label" %in% colnames(JY23@meta.data)) {
  stop("The 'Label' column was not successfully added to the JY23 Seurat object's metadata.")
}

JY23_tumor <- subset(JY23, subset = Label == "Tumor")
JY23_nontumor <- subset(JY23, subset = Label == "Non-Tumor")

# Step 6: Save the tumor and non-tumor Seurat objects for all samples
saveRDS(JY4_tumor, file = "JY4/JY4_tumor.rds")
saveRDS(JY4_nontumor, file = "JY4/JY4_nontumor.rds")

saveRDS(JY8_tumor, file = "JY8/JY8_tumor.rds")
saveRDS(JY8_nontumor, file = "JY8/JY8_nontumor.rds")

saveRDS(JY10_tumor, file = "JY10/JY10_tumor.rds")
saveRDS(JY10_nontumor, file = "JY10/JY10_nontumor.rds")

saveRDS(JY23_tumor, file = "JY23/JY23_tumor.rds")
saveRDS(JY23_nontumor, file = "JY23/JY23_nontumor.rds")

# Step 7: Output summary information for all samples
cat("Summary for JY4:\n")
cat("Number of tumor cells:", ncol(JY4_tumor), "\n")
cat("Number of non-tumor cells:", ncol(JY4_nontumor), "\n\n")

cat("Summary for JY8:\n")
cat("Number of tumor cells:", ncol(JY8_tumor), "\n")
cat("Number of non-tumor cells:", ncol(JY8_nontumor), "\n\n")

cat("Summary for JY10:\n")
cat("Number of tumor cells:", ncol(JY10_tumor), "\n")
cat("Number of non-tumor cells:", ncol(JY10_nontumor), "\n\n")

cat("Summary for JY23:\n")
cat("Number of tumor cells:", ncol(JY23_tumor), "\n")
cat("Number of non-tumor cells:", ncol(JY23_nontumor), "\n")