# =============================================================================
# 01_qc_filtering_integration.R
# -----------------------------------------------------------------------------
# Step 1: from raw 10x counts to an integrated, clustered Seurat object.
#
# DATASETS (the only two used in this chapter):
#   - T21BM_male04 : Trisomy 21 (Down syndrome) fetal bone marrow
#   - D21_male35   : Disomic (typical control)  fetal bone marrow
#
# These are real CellRanger outputs on the shared AWS share (read-only), so you
# run the full Read10X -> QC -> filter -> merge -> integrate -> cluster workflow
# directly. No downloads.
#
# Companion lesson: lessons/01_qc_filtering_integration.md
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(SoupX)       # ambient RNA estimation (demonstrated below)
library(tidyverse)   # dplyr + tidyr + ggplot2, etc.

# -----------------------------------------------------------------------------
# 1. Read the two 10x matrices and make Seurat objects
# -----------------------------------------------------------------------------
# Read10X() expects a CellRanger filtered_feature_bc_matrix/ folder containing
# barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz.
t21.mat <- Read10X(data.dir = file.path(RAW10X_DIR, T21_SAMPLE, "outs", "filtered_feature_bc_matrix"))
d21.mat <- Read10X(data.dir = file.path(RAW10X_DIR, D21_SAMPLE, "outs", "filtered_feature_bc_matrix"))

# project = sample name; this fills orig.ident, which we rely on later.
t21.obj <- CreateSeuratObject(counts = t21.mat, project = "T21BM_male04")
d21.obj <- CreateSeuratObject(counts = d21.mat, project = "D21_male35")

# -----------------------------------------------------------------------------
# 2. (Recommended in real runs) Ambient RNA correction with SoupX
# -----------------------------------------------------------------------------
# SoupX estimates the fraction of counts coming from free-floating ("soup")
# mRNA rather than the cell, then subtracts it. It needs BOTH the raw
# (unfiltered) and filtered matrices plus clusters, which load10X() reads
# straight from the CellRanger outs/ folder. Shown here so you know where it
# fits; left commented so this script runs on the filtered matrices above.
#
# sc  <- load10X(file.path(RAW10X_DIR, T21_SAMPLE, "outs"))  # raw + filtered + clusters
# sc  <- autoEstCont(sc)                                      # estimate contamination
# out <- adjustCounts(sc)                                     # corrected counts
# t21.obj <- CreateSeuratObject(counts = out, project = "T21BM_male04")

# -----------------------------------------------------------------------------
# 3. Add experimental metadata
# -----------------------------------------------------------------------------
t21.obj@meta.data$T21.status <- "T21"
t21.obj@meta.data$gender     <- "male"
d21.obj@meta.data$T21.status <- "D21"
d21.obj@meta.data$gender     <- "male"

# -----------------------------------------------------------------------------
# 4. Quality-control metric: percent mitochondrial reads
# -----------------------------------------------------------------------------
# High mito % usually means a dying/lysed cell whose cytoplasmic mRNA leaked out.
# Human mitochondrial genes start with "MT-" (mouse would be "mt-").
t21.obj[["percent.mt"]] <- PercentageFeatureSet(t21.obj, pattern = "^MT-")
d21.obj[["percent.mt"]] <- PercentageFeatureSet(d21.obj, pattern = "^MT-")

# Inspect the three core QC distributions BEFORE filtering.
VlnPlot(t21.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(d21.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# -----------------------------------------------------------------------------
# 5. Filter out low-quality cells
# -----------------------------------------------------------------------------
# nFeature_RNA > 200   : remove empty droplets / debris (too few genes)
# nFeature_RNA < 2500  : remove likely doublets (suspiciously many genes)
# percent.mt   < 5     : remove dying cells
# NOTE: these thresholds are dataset-dependent. Set them by LOOKING at the
# violin plots above, not by copying numbers blindly.
t21.obj <- subset(t21.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
d21.obj <- subset(d21.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Re-inspect after filtering to confirm the extreme tails are gone.
VlnPlot(t21.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(d21.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# -----------------------------------------------------------------------------
# 6. Merge the two datasets and run the standard preprocessing chain
# -----------------------------------------------------------------------------
merged <- merge(t21.obj, d21.obj)

merged <- NormalizeData(merged)        # log-normalize counts
merged <- FindVariableFeatures(merged) # find the most informative genes
merged <- ScaleData(merged)            # center/scale so PCA isn't dominated by high-expression genes
merged <- RunPCA(merged)               # linear dimensionality reduction

# -----------------------------------------------------------------------------
# 7. Integrate the two samples (remove batch/condition technical differences)
# -----------------------------------------------------------------------------
# RPCA = Reciprocal PCA integration. It aligns shared cell populations across
# the T21 and D21 samples so the same cell type from each sample overlaps,
# instead of forming two separate blobs driven by batch.
obj <- IntegrateLayers(
  object         = merged,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)

# -----------------------------------------------------------------------------
# 8. Cluster and build the UMAP on the INTEGRATED space
# -----------------------------------------------------------------------------
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30,
               reduction.name = "umap.rpca", n.neighbors = 100, min.dist = 0.5)
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "rpca_clusters")

# Visualize clusters and metadata to sanity-check the integration.
DimPlot(obj, pt.size = 1, group.by = "seurat_clusters", reduction = "umap.rpca")
DimPlot(obj, pt.size = 1, group.by = "T21.status",     reduction = "umap.rpca")

# Save to YOUR writable output folder (the share is read-only).
save(obj, file = file.path(OUT_DIR, "integrated_clustered.RData"))
