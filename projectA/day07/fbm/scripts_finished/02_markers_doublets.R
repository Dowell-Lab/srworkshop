# =============================================================================
# 02_markers_doublets.R
# -----------------------------------------------------------------------------
# Step 2: find marker / differentially expressed genes, and detect doublets.
#
# Uses the integrated, clustered object produced by 01_qc_filtering_integration.R
# (built from the T21BM_male04 + D21_male35 datasets only).
#
# Companion lesson: lessons/02_markers_doublets.md
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(patchwork)

load(file.path(OUT_DIR, "integrated_clustered.RData"))  # -> obj

# -----------------------------------------------------------------------------
# 1. Differentially expressed genes ACROSS all clusters
# -----------------------------------------------------------------------------
# JoinLayers collapses the per-sample layers (created during integration) back
# into one matrix so marker testing sees all cells together.
joint.layers.obj <- JoinLayers(obj)

# FindAllMarkers: for each cluster, genes up in that cluster vs all others.
all.clusters.diff.genes <- FindAllMarkers(joint.layers.obj, assay = "RNA")

# Top 10 markers per cluster by fold change (the genes that DEFINE each cluster).
all.clusters.diff.genes %>%
  group_by(cluster) %>%
  top_n(wt = avg_log2FC, n = 10)

# Top 10 by pct.1 = fraction of cells in the cluster expressing the gene
# (useful for finding clean, highly-detected markers).
all.clusters.diff.genes %>%
  group_by(cluster) %>%
  top_n(wt = pct.1, n = 10)

# -----------------------------------------------------------------------------
# 2. Differentially expressed genes BETWEEN two specific clusters
# -----------------------------------------------------------------------------
# ident.1 vs ident.2: a direct pairwise comparison (positive log2FC = up in ident.1).
cluster1.cluster2.diff.genes <- FindMarkers(joint.layers.obj, ident.1 = "1", ident.2 = "2")
head(cluster1.cluster2.diff.genes)

# -----------------------------------------------------------------------------
# 3. Differentially expressed genes BY METADATA (condition), not cluster
# -----------------------------------------------------------------------------
# Switch the "identity" of each cell to its T21 status, then test T21 vs D21.
joint.layers.obj <- SetIdent(joint.layers.obj, value = "T21.status")
diff.genes.by.status <- FindAllMarkers(joint.layers.obj)

# Condition comparison WITHIN a single cluster: subset to cluster 1, then test
# T21 vs D21 inside it (controls for cell-type composition differences).
joint.layers.obj <- SetIdent(joint.layers.obj, value = "seurat_clusters")
cluster1.subset   <- subset(joint.layers.obj, idents = "1")
cluster1.subset   <- SetIdent(cluster1.subset, value = "T21.status")
diff.genes.cluster1.status <- FindAllMarkers(cluster1.subset)
diff.genes.cluster1.status

# -----------------------------------------------------------------------------
# 4. Visualize where a marker gene is expressed
# -----------------------------------------------------------------------------
DefaultAssay(joint.layers.obj) <- "RNA"
FeaturePlot(joint.layers.obj, features = "SLC4A1",  reduction = "umap.rpca")  # erythroid marker
FeaturePlot(joint.layers.obj, features = "FCGR3B",  reduction = "umap.rpca")  # neutrophil marker

# Save the marker tables for downstream use.
save(all.clusters.diff.genes, file = file.path(OUT_DIR, "all_clusters_diff_genes.RData"))

# -----------------------------------------------------------------------------
# 5. Doublet detection with DoubletFinder
# -----------------------------------------------------------------------------
# A "doublet" is two cells captured in one droplet, producing a fake hybrid
# profile. DoubletFinder works by spiking in artificial doublets and asking
# which real cells look most like them.
#   pN   = proportion of artificial doublets to generate
#   pK   = neighborhood size (ideally tuned with paramSweep; 0.09 used here)
#   nExp = how many doublets we EXPECT (here ~10% of cells)
#   PCs  = principal components to use
nExp <- round(ncol(joint.layers.obj) * 0.10)

joint.layers.obj <- doubletFinder(
  joint.layers.obj,
  pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = FALSE
)

# The classification lands in a metadata column named "DF.classifications_...".
DF.name <- colnames(joint.layers.obj@meta.data)[grepl("DF.classification", colnames(joint.layers.obj@meta.data))]

# Color the UMAP by singlet vs doublet.
DimPlot(joint.layers.obj, group.by = DF.name, reduction = "umap.rpca")

# Typically you would now drop the doublets:
# singlets <- subset(joint.layers.obj, subset = !!sym(DF.name) == "Singlet")
