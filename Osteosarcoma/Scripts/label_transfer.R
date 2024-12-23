library(Seurat)
library(tidyverse)

#' Process and validate input data
#' @param reference_obj Reference Seurat object
#' @param query_obj Query Seurat object
#' @param annotations Annotation data frame
#' @param assay_name Assay name to use
#' @return List of processed objects
validate_inputs <- function(reference_obj, query_obj, annotations, assay_name) {
  if (!all(c(assay_name) %in% Assays(reference_obj)) || 
      !all(c(assay_name) %in% Assays(query_obj))) {
    stop("Required assays not found")
  }
  
  common_barcodes <- intersect(colnames(reference_obj), rownames(annotations))
  if (length(common_barcodes) == 0) {
    stop("No matching barcodes between reference and annotations")
  }
  
  reference_obj$Graph.based <- NA
  reference_obj$Graph.based[common_barcodes] <- annotations[common_barcodes, "Graph.based"]
  reference_subset <- subset(reference_obj, subset = !is.na(Graph.based) & Graph.based != "")
  
  list(reference = reference_subset, query = query_obj)
}

#' Process data for transfer
#' @param object Seurat object
#' @param variable_features Number of variable features
#' @param assay_name Assay name
#' @return Processed Seurat object
process_data <- function(object, variable_features, assay_name) {
  DefaultAssay(object) <- assay_name
  object <- NormalizeData(object)
  if ("reference" %in% attributes(object)$names) {
    object <- FindVariableFeatures(object, nfeatures = variable_features)
  }
  object <- ScaleData(object, features = VariableFeatures(object))
  if ("reference" %in% attributes(object)$names) {
    object <- RunPCA(object, features = VariableFeatures(object))
  }
  object
}

#' Perform label transfer
#' @param reference Processed reference object
#' @param query Processed query object
#' @param params Transfer parameters
#' @return Query object with predictions
transfer_labels <- function(reference, query, params) {
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    reduction = "pcaproject",
    features = VariableFeatures(reference),
    reference.assay = params$assay_name,
    query.assay = params$assay_name,
    normalization.method = "LogNormalize"
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference$Graph.based,
    dims = params$pca_dims
  )
  
  query$predicted.celltype <- predictions$predicted.id
  query$prediction.score <- predictions$prediction.score.max
  query$predicted.celltype.filtered <- ifelse(
    query$prediction.score >= params$confidence_threshold,
    query$predicted.celltype,
    "Unlabeled"
  )
  
  query
}

#' Create non-OS annotations
#' @param annotations Original annotations
#' @return Modified annotations
create_non_os_annotations <- function(annotations) {
  non_os <- annotations
  non_os$Graph.based[non_os$Graph.based == "Osteoblastic_OS_cells"] <- ""
  non_os
}

#' Combine predictions from both transfers
#' @param non_os_query Query with non-OS predictions
#' @param full_query Query with full predictions
#' @param threshold Confidence threshold
#' @return List with combined results
combine_predictions <- function(non_os_query, full_query, threshold) {
  combined_meta <- data.frame(
    spot_id = colnames(non_os_query),
    non_os_prediction = non_os_query$predicted.celltype.filtered,
    non_os_confidence = non_os_query$prediction.score,
    full_prediction = full_query$predicted.celltype,
    full_confidence = full_query$prediction.score,
    row.names = colnames(non_os_query)
  )
  
  combined_meta$final_prediction <- combined_meta$non_os_prediction
  unlabeled_mask <- combined_meta$non_os_prediction == "Unlabeled"
  combined_meta$final_prediction[unlabeled_mask] <- ifelse(
    combined_meta$full_confidence[unlabeled_mask] >= threshold,
    combined_meta$full_prediction[unlabeled_mask],
    "Unlabeled"
  )
  
  combined_meta$final_confidence <- combined_meta$non_os_confidence
  combined_meta$final_confidence[unlabeled_mask] <- 
    combined_meta$full_confidence[unlabeled_mask]
  
  non_os_query$final.celltype <- combined_meta$final_prediction
  non_os_query$final.confidence <- combined_meta$final_confidence
  
  list(
    object = non_os_query,
    metadata = combined_meta
  )
}

#' Main pipeline function
#' @param reference_path Path to reference RDS
#' @param query_path Path to query RDS
#' @param annotations_path Path to annotations CSV
#' @param output_dir Output directory
#' @param params Analysis parameters
#' @return Combined results list
run_label_transfer <- function(
    reference_path,
    query_path,
    annotations_path,
    output_dir,
    params = list(
      confidence_threshold = 0.7,
      variable_features = 2000,
      pca_dims = 1:30,
      assay_name = "Spatial.008um"
    )
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  reference <- readRDS(reference_path)
  query <- readRDS(query_path)
  annotations <- read.csv(annotations_path, row.names = 1)
  
  # Create and validate data versions
  non_os_annotations <- create_non_os_annotations(annotations)
  non_os_data <- validate_inputs(reference, query, non_os_annotations, params$assay_name)
  full_data <- validate_inputs(reference, query, annotations, params$assay_name)
  
  # Process data
  non_os_reference <- process_data(non_os_data$reference, params$variable_features, params$assay_name)
  full_reference <- process_data(full_data$reference, params$variable_features, params$assay_name)
  query_processed <- process_data(query, params$variable_features, params$assay_name)
  
  # Perform transfers
  non_os_results <- transfer_labels(non_os_reference, query_processed, params)
  full_results <- transfer_labels(full_reference, query_processed, params)
  
  # Combine results
  results <- combine_predictions(non_os_results, full_results, params$confidence_threshold)
  
  # Save outputs
  saveRDS(results$object, file.path(output_dir, "annotated_query.rds"))
  write.csv(results$metadata, file.path(output_dir, "predictions.csv"))
  
  results
}

# Example usage:
# results <- run_label_transfer(
#   reference_path = "rds_files/JY10_object.rds",
#   query_path = "rds_files/JY8_object.rds",
#   annotations_path = "barcodes/JY10_ANNOTATED_ORIGINAL.csv",
#   output_dir = "output/JY8_transfer"
# )