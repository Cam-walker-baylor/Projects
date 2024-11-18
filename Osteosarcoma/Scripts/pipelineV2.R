library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

marker_list <- list(
  "Osteoblastic OS cells" = c("COL1A1", "CDH11", "RUNX2"),
  "Proliferating osteoblastic OS cells" = c("TOP2A", "PCNA", "MKI67"),
  "Chondroblastic OS cells" = c("ACAN", "COL2A1", "SOX9"),
  "Osteoclastic cells" = c("CTSK", "MMP9"),
  "TILs (T and NK cells)" = c("IL7R", "CD3D", "NKG7"),
  "Myeloid cells" = c("CD74", "CD14", "FCGR3A"),
  "Fibroblasts" = c("COL1A1", "LUM", "DCN"),
  "Pericytes" = c("ACTA2", "RGS5"),
  "MSCs" = c("CXCL12", "SFRP2", "MME"),
  "Myoblasts" = c("MYLPF", "MYL1"),
  "Endothelial cells" = c("PECAM1", "VWF")
)


JY4_tumor <- readRDS(file = "JY4_tumor.rds")
JY4_nontumor <- readRDS(file = "JY4_nontumor.rds")
JY8_tumor <- readRDS(file = "JY8_tumor.rds")
JY8_nontumor <- readRDS(file = "JY8_nontumor.rds")
JY10_tumor <- readRDS(file = "JY10_tumor.rds")
JY10_nontumor <- readRDS(file = "JY10_nontumor.rds")
JY23_tumor <- readRDS(file = "JY23_tumor.rds")
JY23_nontumor <- readRDS(file = "JY23_nontumor.rds")


# normalize both 8um and 16um bins
DefaultAssay(JY4_tumor) <- "Spatial.008um"
JY4_tumor <- NormalizeData(JY4_tumor)
DefaultAssay(JY4_nontumor) <- "Spatial.008um"
JY4_nontumor <- NormalizeData(JY4_nontumor)

DefaultAssay(JY8_tumor) <- "Spatial.008um"
JY8_tumor <- NormalizeData(JY8_tumor)
DefaultAssay(JY8_nontumor) <- "Spatial.008um"
JY8_nontumor <- NormalizeData(JY8_nontumor)

DefaultAssay(JY10_tumor) <- "Spatial.008um"
JY10_tumor <- NormalizeData(JY10_tumor)
DefaultAssay(JY10_nontumor) <- "Spatial.008um"
JY10_nontumor <- NormalizeData(JY10_nontumor)

DefaultAssay(JY23_tumor) <- "Spatial.008um"
JY23_tumor <- NormalizeData(JY23_tumor)
DefaultAssay(JY23_nontumor) <- "Spatial.008um"
JY23_nontumor <- NormalizeData(JY23_nontumor)


# note that data is already normalized
DefaultAssay(JY4_tumor) <- "Spatial.008um"
JY4_tumor <- FindVariableFeatures(JY4_tumor)
JY4_tumor <- ScaleData(JY4_tumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY4_tumor <- SketchData(
  object = JY4_tumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY4_nontumor) <- "Spatial.008um"
JY4_nontumor <- FindVariableFeatures(JY4_nontumor)
JY4_nontumor <- ScaleData(JY4_nontumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY4_nontumor <- SketchData(
  object = JY4_nontumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY8_tumor) <- "Spatial.008um"
JY8_tumor <- FindVariableFeatures(JY8_tumor)
JY8_tumor <- ScaleData(JY8_tumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY8_tumor <- SketchData(
  object = JY8_tumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY8_nontumor) <- "Spatial.008um"
JY8_nontumor <- FindVariableFeatures(JY8_nontumor)
JY8_nontumor <- ScaleData(JY8_nontumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY8_nontumor <- SketchData(
  object = JY8_nontumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY10_tumor) <- "Spatial.008um"
JY10_tumor <- FindVariableFeatures(JY10_tumor)
JY10_tumor <- ScaleData(JY10_tumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY10_tumor <- SketchData(
  object = JY10_tumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY10_nontumor) <- "Spatial.008um"
JY10_nontumor <- FindVariableFeatures(JY10_nontumor)
JY10_nontumor <- ScaleData(JY10_nontumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY10_nontumor <- SketchData(
  object = JY10_nontumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY23_tumor) <- "Spatial.008um"
JY23_tumor <- FindVariableFeatures(JY23_tumor)
JY23_tumor <- ScaleData(JY23_tumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY23_tumor <- SketchData(
  object = JY23_tumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(JY23_nontumor) <- "Spatial.008um"
JY23_nontumor <- FindVariableFeatures(JY23_nontumor)
JY23_nontumor <- ScaleData(JY23_nontumor)
# we select 50,0000 cells and create a new 'sketch' assay
JY23_nontumor <- SketchData(
  object = JY23_nontumor, #object stays object
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(JY4_tumor) <- "sketch"
DefaultAssay(JY4_nontumor) <- "sketch"

DefaultAssay(JY8_tumor) <- "sketch"
DefaultAssay(JY8_nontumor) <- "sketch"

DefaultAssay(JY10_tumor) <- "sketch"
DefaultAssay(JY10_nontumor) <- "sketch"

DefaultAssay(JY23_tumor) <- "sketch"
DefaultAssay(JY23_nontumor) <- "sketch"


# perform clustering workflow
JY4_tumor <- FindVariableFeatures(JY4_tumor)
JY4_tumor <- ScaleData(JY4_tumor)
JY4_tumor <- RunPCA(JY4_tumor, assay = "sketch", reduction.name = "pca.sketch")

JY4_nontumor <- FindVariableFeatures(JY4_nontumor)
JY4_nontumor <- ScaleData(JY4_nontumor)
JY4_nontumor <- RunPCA(JY4_nontumor, assay = "sketch", reduction.name = "pca.sketch")

JY8_tumor <- FindVariableFeatures(JY8_tumor)
JY8_tumor <- ScaleData(JY8_tumor)
JY8_tumor <- RunPCA(JY8_tumor, assay = "sketch", reduction.name = "pca.sketch")

JY8_nontumor <- FindVariableFeatures(JY8_nontumor)
JY8_nontumor <- ScaleData(JY8_nontumor)
JY8_nontumor <- RunPCA(JY8_nontumor, assay = "sketch", reduction.name = "pca.sketch")

JY10_tumor <- FindVariableFeatures(JY10_tumor)
JY10_tumor <- ScaleData(JY10_tumor)
JY10_tumor <- RunPCA(JY10_tumor, assay = "sketch", reduction.name = "pca.sketch")

JY10_nontumor <- FindVariableFeatures(JY10_nontumor)
JY10_nontumor <- ScaleData(JY10_nontumor)
JY10_nontumor <- RunPCA(JY10_nontumor, assay = "sketch", reduction.name = "pca.sketch")

JY23_tumor <- FindVariableFeatures(JY23_tumor)
JY23_tumor <- ScaleData(JY23_tumor)
JY23_tumor <- RunPCA(JY23_tumor, assay = "sketch", reduction.name = "pca.sketch")

JY23_nontumor <- FindVariableFeatures(JY23_nontumor)
JY23_nontumor <- ScaleData(JY23_nontumor)
JY23_nontumor <- RunPCA(JY23_nontumor, assay = "sketch", reduction.name = "pca.sketch")

# Determine the 'dimensionality' of the dataset
ElbowPlot(JY4_tumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY4_tumor.png", plot = ElbowPlot(JY4_tumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY4_nontumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY4_nontumor.png", plot = ElbowPlot(JY4_nontumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY8_tumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY8_tumor.png", plot = ElbowPlot(JY8_tumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY8_nontumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY8_nontumor.png", plot = ElbowPlot(JY8_nontumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY10_tumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY10_tumor.png", plot = ElbowPlot(JY10_tumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY10_nontumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY10_nontumor.png", plot = ElbowPlot(JY10_nontumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY23_tumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY23_tumor.png", plot = ElbowPlot(JY23_tumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)

ElbowPlot(JY23_nontumor, reduction = "pca.sketch")
ggsave(filename = "ElbowPlot_JY23_nontumor.png", plot = ElbowPlot(JY23_nontumor, reduction = "pca.sketch"), width = 8, height = 6, dpi = 300)


JY4_tumor <- FindNeighbors(JY4_tumor, assay = "sketch", reduction = "pca.sketch", dims = 1:7)
JY4_tumor <- FindClusters(JY4_tumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY4_tumor <- RunUMAP(JY4_tumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:7)

JY4_nontumor <- FindNeighbors(JY4_nontumor, assay = "sketch", reduction = "pca.sketch", dims = 1:4)
JY4_nontumor <- FindClusters(JY4_nontumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY4_nontumor <- RunUMAP(JY4_nontumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:4)

JY8_tumor <- FindNeighbors(JY8_tumor, assay = "sketch", reduction = "pca.sketch", dims = 1:7)
JY8_tumor <- FindClusters(JY8_tumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY8_tumor <- RunUMAP(JY8_tumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:7)

JY8_nontumor <- FindNeighbors(JY8_nontumor, assay = "sketch", reduction = "pca.sketch", dims = 1:4)
JY8_nontumor <- FindClusters(JY8_nontumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY8_nontumor <- RunUMAP(JY8_nontumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:4)

JY10_tumor <- FindNeighbors(JY10_tumor, assay = "sketch", reduction = "pca.sketch", dims = 1:7)
JY10_tumor <- FindClusters(JY10_tumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY10_tumor <- RunUMAP(JY10_tumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:7)

JY10_nontumor <- FindNeighbors(JY10_nontumor, assay = "sketch", reduction = "pca.sketch", dims = 1:4)
JY10_nontumor <- FindClusters(JY10_nontumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY10_nontumor <- RunUMAP(JY10_nontumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:4)

JY23_tumor <- FindNeighbors(JY23_tumor, assay = "sketch", reduction = "pca.sketch", dims = 1:7)
JY23_tumor <- FindClusters(JY23_tumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY23_tumor <- RunUMAP(JY23_tumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:7)

JY23_nontumor <- FindNeighbors(JY23_nontumor, assay = "sketch", reduction = "pca.sketch", dims = 1:4)
JY23_nontumor <- FindClusters(JY23_nontumor, cluster.name = "seurat_cluster.sketched", resolution = 1)
JY23_nontumor <- RunUMAP(JY23_nontumor, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:4)

# Filter out cells without a cluster assignment in seurat_cluster.sketched
JY4_tumor <- subset(JY4_tumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")
JY4_nontumor <- subset(JY4_nontumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")

JY8_tumor <- subset(JY8_tumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")
JY8_nontumor <- subset(JY8_nontumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")

JY10_tumor <- subset(JY10_tumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")
JY10_nontumor <- subset(JY10_nontumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")

JY23_tumor <- subset(JY23_tumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")
JY23_nontumor <- subset(JY23_nontumor, subset = !is.na(seurat_cluster.sketched) & seurat_cluster.sketched != "")


JY4_tumor <- ProjectData(
  object = JY4_tumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:3,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY4_nontumor <- ProjectData(
  object = JY4_nontumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:2,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY8_tumor <- ProjectData(
  object = JY8_tumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:3,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY8_nontumor <- ProjectData(
  object = JY8_nontumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:2,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY10_tumor <- ProjectData(
  object = JY10_tumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:3,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY10_nontumor <- ProjectData(
  object = JY10_nontumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:2,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY23_tumor <- ProjectData(
  object = JY23_tumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:3,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

JY23_nontumor <- ProjectData(
  object = JY23_nontumor,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:2,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# Visualize the clustering results
JY4_tumor_umap_plot <- DimPlot(JY4_tumor, reduction = "umap.sketch")
JY4_nontumor_umap_plot <- DimPlot(JY4_nontumor, reduction = "umap.sketch")

JY8_tumor_umap_plot <- DimPlot(JY8_tumor, reduction = "umap.sketch")
JY8_nontumor_umap_plot <- DimPlot(JY8_nontumor, reduction = "umap.sketch")

JY10_tumor_umap_plot <- DimPlot(JY10_tumor, reduction = "umap.sketch")
JY10_nontumor_umap_plot <- DimPlot(JY10_nontumor, reduction = "umap.sketch")

JY23_tumor_umap_plot <- DimPlot(JY23_tumor, reduction = "umap.sketch")
JY23_nontumor_umap_plot <- DimPlot(JY23_nontumor, reduction = "umap.sketch")


# Save the plot
ggsave(filename = "DimPlot_JY4_tumor.png", plot = JY4_tumor_umap_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = "DimPlot_JY4_nontumor.png", plot = JY4_nontumor_umap_plot, width = 8, height = 6, dpi = 300)

ggsave(filename = "DimPlot_JY8_tumor.png", plot = JY8_tumor_umap_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = "DimPlot_JY8_nontumor.png", plot = JY8_nontumor_umap_plot, width = 8, height = 6, dpi = 300)

ggsave(filename = "DimPlot_JY10_tumor.png", plot = JY10_tumor_umap_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = "DimPlot_JY10_nontumor.png", plot = JY10_nontumor_umap_plot, width = 8, height = 6, dpi = 300)

ggsave(filename = "DimPlot_JY23_tumor.png", plot = JY23_tumor_umap_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = "DimPlot_JY23_nontumor.png", plot = JY23_nontumor_umap_plot, width = 8, height = 6, dpi = 300)


# (Optional) Visualize the spatial plot again with adjustments
JY4_tumor_spatial_plot <- SpatialDimPlot(JY4_tumor, 
                                   images = "slice1.008um", 
                                   alpha = 0.5,           # Adjust transparency
                                   pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY4_tumor_spatial_plot_008.png", plot = JY4_tumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY4_nontumor_spatial_plot <- SpatialDimPlot(JY4_nontumor, 
                                   images = "slice1.008um", 
                                   alpha = 0.5,           # Adjust transparency
                                   pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY4_nontumor_spatial_plot_008.png", plot = JY4_nontumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY8_tumor_spatial_plot <- SpatialDimPlot(JY8_tumor, 
                                          images = "slice1.008um", 
                                          alpha = 0.5,           # Adjust transparency
                                          pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY8_tumor_spatial_plot_008.png", plot = JY8_tumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY8_nontumor_spatial_plot <- SpatialDimPlot(JY8_nontumor, 
                                             images = "slice1.008um", 
                                             alpha = 0.5,           # Adjust transparency
                                             pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY8_nontumor_spatial_plot_008.png", plot = JY8_nontumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY10_tumor_spatial_plot <- SpatialDimPlot(JY10_tumor, 
                                          images = "slice1.008um", 
                                          alpha = 0.5,           # Adjust transparency
                                          pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY10_tumor_spatial_plot_008.png", plot = JY10_tumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY10_nontumor_spatial_plot <- SpatialDimPlot(JY10_nontumor, 
                                             images = "slice1.008um", 
                                             alpha = 0.5,           # Adjust transparency
                                             pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY10_nontumor_spatial_plot_008.png", plot = JY10_nontumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY23_tumor_spatial_plot <- SpatialDimPlot(JY23_tumor, 
                                          images = "slice1.008um", 
                                          alpha = 0.5,           # Adjust transparency
                                          pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY23_tumor_spatial_plot_008.png", plot = JY23_tumor_spatial_plot, width = 8, height = 6, dpi = 300)

JY23_nontumor_spatial_plot <- SpatialDimPlot(JY23_nontumor, 
                                             images = "slice1.008um", 
                                             alpha = 0.5,           # Adjust transparency
                                             pt.size.factor = 15.0   # Adjust spot size
)

ggsave("JY23_nontumor_spatial_plot_008.png", plot = JY23_nontumor_spatial_plot, width = 8, height = 6, dpi = 300)

# Optional only needed if cells need to be filtered out
## JY4_object_filtered <- subset(JY4_object, idents = WhichCells(JY4_object, idents = which(table(Idents(JY4_object)) > 100)))

# Find markers for each cluster
JY4_tumor_cluster_markers <- FindAllMarkers(JY4_tumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
JY4_nontumor_cluster_markers <- FindAllMarkers(JY4_nontumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

JY8_tumor_cluster_markers <- FindAllMarkers(JY8_tumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
JY8_nontumor_cluster_markers <- FindAllMarkers(JY8_nontumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

JY10_tumor_cluster_markers <- FindAllMarkers(JY10_tumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
JY10_nontumor_cluster_markers <- FindAllMarkers(JY10_nontumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

JY23_tumor_cluster_markers <- FindAllMarkers(JY23_tumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
JY23_nontumor_cluster_markers <- FindAllMarkers(JY23_nontumor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

# For tumor clusters
tumor_markers <- JY10_tumor_cluster_markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(tumor_markers, "tumor_cluster_markers.csv", row.names = FALSE)

# For non-tumor clusters  
nontumor_markers <- JY10_nontumor_cluster_markers %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(nontumor_markers, "nontumor_cluster_markers.csv", row.names = FALSE)

####################################


# 9. Annotate cell types using marker list
# Assuming marker_list is a named list of markers for each cell type
JY4_annotate_clusters <- function(seurat_object, JY4_cluster_markers, marker_list) {
  annotations <- list()
  for (cluster in levels(seurat_object@active.ident)) {
    cluster_genes <- JY4_cluster_markers %>%
      filter(cluster == !!cluster) %>%
      top_n(n = 1000, wt = avg_log2FC) %>%
      pull(gene)
    
    scores <- sapply(marker_list, function(markers) {
      sum(markers %in% cluster_genes) / length(markers)
    })
    
    # Assign the highest scoring cell type or "Unlabeled" if no match is found
    if (max(scores) > 0) {
      annotations[[cluster]] <- names(which.max(scores))
    } else {
      annotations[[cluster]] <- "Unlabeled"  # Assign "Unlabeled" if no match is found
    }
  }
  
  new_idents <- unlist(annotations)[as.character(seurat_object@active.ident)]
  seurat_object <- RenameIdents(seurat_object, new_idents)
  return(seurat_object)
  
}

JY10_annotate_clusters <- function(seurat_object, JY10_cluster_markers, marker_list) {
  annotations <- list()
  for (cluster in levels(seurat_object@active.ident)) {
    cluster_genes <- JY10_cluster_markers %>%
      filter(cluster == !!cluster) %>%
      top_n(n = 25, wt = avg_log2FC) %>%
      pull(gene)
    
    scores <- sapply(marker_list, function(markers) {
      sum(markers %in% cluster_genes) / length(markers)
    })
    
    # Assign the highest scoring cell type or "Unlabeled" if no match is found
    if (max(scores) > 0) {
      annotations[[cluster]] <- names(which.max(scores))
    } else {
      annotations[[cluster]] <- "Unlabeled"  # Assign "Unlabeled" if no match is found
    }
  }
  
  new_idents <- unlist(annotations)[as.character(seurat_object@active.ident)]
  seurat_object <- RenameIdents(seurat_object, new_idents)
  return(seurat_object)
  
}


# Apply the annotation function
JY4_object <- JY4_annotate_clusters(JY4_object, JY4_cluster_markers, marker_list)
JY10_object <- JY10_annotate_clusters(JY10_object, JY10_cluster_markers, marker_list)


# Check the identities after annotation
table(Idents(JY4_object))
table(Idents(JY10_object))

# 10. Visualize the annotated clusters
JY4_annotated_plot <- DimPlot(JY4_object, reduction = "umap.sketch", label = TRUE, pt.size = 0.5) + NoLegend()

ggsave(filename = "JY4_annotated_plot.png", plot = JY4_annotated_plot, width = 8, height = 6, dpi = 300)


# 11. Save the results
saveRDS(JY4_object, file = "annotated_visium_hd_object.rds")


# Filter marker genes based on what is present in JY4_object
JY4_filtered_marker_list <- lapply(marker_list, function(genes) {
  genes[genes %in% rownames(JY4_object_subset)]
})

JY10_filtered_marker_list <- lapply(marker_list, function(genes) {
  genes[genes %in% rownames(JY10_object_subset)]
})
