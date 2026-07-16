# =============================================================================
# GATA1 track - Script 03: Multi-reference SingleR annotation + marker genes
# -----------------------------------------------------------------------------
# Assigns a cell-type label to every cell automatically with SingleR, and
# shows WHY you might consult more than one reference atlas for hematopoiesis.
# Then sanity-checks the labels against canonical marker genes.
#
# Inputs : OUT_DIR/gata1_combined_clustered.rds  (from script 02)
# Outputs: an annotated object + UMAPs / FeaturePlots in OUT_DIR.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(patchwork)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_clustered_subsampled.rds"))
#If you didn't make this yet, use mine!
#combined <- readRDS(file.path(COOKING, "gata1_combined_clustered_subsampled.rds"))


# ---- 1. Prepare a single expression matrix for SingleR ----------------------
# After merge(), the RNA assay holds one "data" layer PER sample
# (data.EuploidGATA1sD7, data.T21wtGATA1D9, ...). SingleR wants ONE matrix,
# so JoinLayers stitches them into a unified "data" layer first.
Layers(combined[["RNA"]])           # inspect the per-sample layers
combined_joined <- JoinLayers(combined)


sce_counts <- LayerData(combined_joined, assay = "RNA", layer = "data")  # log-normalized


# ---- 2. Choose reference atlases --------------------------------------------
# celldex ships several curated references. Which one fits depends on the
# tissue. For iPSC-derived blood you have a few sensible options:
#
#
#
#   BlueprintEncodeData()            - broad human stroma + immune
#   HumanPrimaryCellAtlasData()      - broad human primary cells (our default)
#   MonacoImmuneData()               - fine-grained immune subsets
#   NovershternHematopoieticData()   - classic hematopoietic hierarchy
#   DatabaseImmuneCellExpressionData() - DICE immune populations
#
# searchReferences() lists more. We annotate with the broad atlas first, then
# (optionally) a blood-specific one, and compare. Loading a reference pulls it
# from the celldex cache; it is not part of the read-only data share.

#!!!!!!!!!!          make sure to finish with NovershternHematopoieticData()    !!!!!

#whichlabelref <- "BlueprintEncodeData"
#ref_broad <- BlueprintEncodeData()

#whichlabelref <- "HumanPrimaryCellAtlasData"
#ref_broad <- HumanPrimaryCellAtlasData()

#whichlabelref <- "MonacoImmuneData"
#ref_broad <- MonacoImmuneData()

#whichlabelref <- "DatabaseImmuneCellExpressionData"
#ref_broad <- DatabaseImmuneCellExpressionData()

whichlabelref <- "NovershternHematopoieticData"
ref_broad <- NovershternHematopoieticData()



#If you have many CPU
#run_singler <- function(ref, test_mat) {
#  SingleR(
#    test   = test_mat,
#    ref    = assay(ref, "logcounts"),
#    labels = ref$label.main
#  )
#}

#If you don't have many CPU

run_singler <- function(ref, test_mat) {
  SingleR(
    test            = test_mat,
    ref             = ref,
    labels          = ref$label.main,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts",
    BPPARAM         = BiocParallel::SerialParam()
  )
}

pred_broad <- run_singler(ref_broad, sce_counts)
table(pred_broad$labels)

# ---- 3. Collapse rare labels into "other" -----------------------------------
# With many references you get a long tail of labels assigned to only a handful
# of cells. Those clutter the UMAP. We keep labels with > 150 cells and bucket
# the rest as "other" for a readable plot, while keeping the full label too.
pred_df <- as.data.frame(pred_broad) %>%
  add_count(labels, name = "label_n") %>%
  mutate(labels_other = case_when(label_n > 150 ~ labels, TRUE ~ "other"))

combined$SingleR_label        <- pred_df$labels
combined$SingleR_labels_other <- pred_df$labels_other
combined$SingleR_pruned       <- pred_df$pruned.labels  # NA where the call was weak
combined$SingleR_delta        <- pred_df$delta.next     # gap to the runner-up label
combined_joined$SingleR_labels_other <- pred_df$labels_other 
# delta.next is a confidence proxy: a small gap means the top two labels were
# nearly tied, so treat that cell's label with caution.

# ---- 4. (Optional) second reference for comparison --------------------------
# Uncomment to annotate with a blood-specific atlas and compare agreement.
# ref_blood  <- NovershternHematopoieticData()
# pred_blood <- run_singler(ref_blood, sce_counts)
# combined$SingleR_label_blood <- pred_blood$labels
# print(table(broad = combined$SingleR_label, blood = combined$SingleR_label_blood))




#Differentially expressed genes BETWEEN two specific clusters

# ident.1 vs ident.2: a direct pairwise comparison (positive log2FC = up in ident.1).
#cluster1.cluster2.diff.genes <- FindMarkers(joint.layers.obj, ident.1 = "1", ident.2 = "2")
#head(cluster1.cluster2.diff.genes)

# -----------------------------------------------------------------------------
# 3. Differentially expressed genes BY METADATA (condition), not cluster
# -----------------------------------------------------------------------------
# Switch the "identity" of each cell to its T21 status, then test T21 vs D21.
#joint.layers.obj <- SetIdent(joint.layers.obj, value = "T21.status")
#diff.genes.by.status <- FindAllMarkers(joint.layers.obj)

# Condition comparison WITHIN a single cluster: subset to cluster 1, then test
# T21 vs D21 inside it (controls for cell-type composition differences).
#joint.layers.obj <- SetIdent(joint.layers.obj, value = "seurat_clusters")
#cluster1.subset   <- subset(joint.layers.obj, idents = "1")
#cluster1.subset   <- SetIdent(cluster1.subset, value = "T21.status")
#diff.genes.cluster1.status <- FindAllMarkers(cluster1.subset)
#diff.genes.cluster1.status



# ---- 5. Visualize the automated annotation ----------------------------------
save_dim <- function(p, file, w = 8, h = 6) {
  ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
}

p<- DimPlot(combined, reduction = "umap", group.by ="SingleR_label", label = TRUE)
p

fn = paste0(whichlabelref, "gata1_umap_singler_label.png")
save_dim(p,
         fn)



p<- DimPlot(combined, reduction = "umap", group.by = "SingleR_pruned", label = TRUE)
p
fn <- paste0(whichlabelref, "gata1_umap_singler_pruned.png")

save_dim(p,
         fn)
p<- DimPlot(combined, reduction = "umap", group.by = "SingleR_labels_other", label = TRUE)
p
fn <- paste0(whichlabelref, "gata1_umap_singler_other.png")
save_dim(p,
         fn)

# ---- 6. Sanity-check with canonical markers ---------------------------------
# Automated labels are a hypothesis. Confirm them with genes you trust:
#   HBG2, HBG1 - fetal hemoglobin (erythroid)
#   GATA1      - master erythroid/megakaryocyte TF (the gene of interest!)
#   GFI1B      - erythroid/megakaryocyte
#   TFRC (CD71)- erythroid progenitors
#   CD34       - hematopoietic stem/progenitor
#   PTPRC (CD45)- pan-leukocyte
#   ITGA4, VAMP8, DIAPH3 - additional lineage markers

agene = "HBG2"

p_markers <- FeaturePlot(combined, features = agene, reduction = "umap",
                         pt.size = 0.3, ncol = 3)
p_markers
fn = paste0("gata1dataset_", agene,".png")
save_dim(p_markers, fn, w = 12, h = 14)


markers <- c("HBG2", "HBG1", "GATA1", "GFI1B", "TFRC",
             "CD34", "PTPRC", "ITGA4", "VAMP8", "DIAPH3")

p_markers <- FeaturePlot(combined, features = markers, reduction = "umap",
                         pt.size = 0.3, ncol = 3)
p_markers
save_dim(p_markers, "gata1_feature_markers.png", w = 12, h = 14)

# Side-by-side: labels next to a key marker.
p_combo <- DimPlot(combined, reduction = "umap", group.by = "SingleR_label", label = TRUE) +
  FeaturePlot(combined, features = "HBG1", reduction = "umap")
p_combo
save_dim(p_combo, "gata1_labels_vs_HBG1.png", w = 13, h = 6)

# ---- 7. Save the fully annotated object -------------------------------------
#saveRDS(combined, file.path(OUT_DIR, "gata1_combined_annotated.rds"))
#message("Saved: ", file.path(OUT_DIR, "gata1_combined_annotated.rds"))

#saveRDS(combined_joined, file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
#message("Saved: ", file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))

# ---- 8. (Optional) What if I don't like what they called? --------------------------

# FindAllMarkers: for each cluster, genes up in that cluster vs all others.
#plan(sequential)

p <- DimPlot(
  combined_joined,
  reduction = "umap",
  group.by = c("SingleR_labels_other", "seurat_clusters"),
  label = TRUE
)
p

Idents(combined_joined)

#or just find some markers

clusters_of_interest <- c(5, 7)  # or character, depending on Idents()


markers_list <- lapply(
  clusters_of_interest,
  function(cl) {
    res <- FindMarkers(
      combined_joined,
      ident.1 = cl,
      assay   = "RNA",
      logfc.threshold = 0.5,
      min.pct        = 0.1
    )
    res$cluster <- cl
    res$gene    <- rownames(res)
    res
  }
)

selected_markers <- do.call(rbind, markers_list)

head(selected_markers)

markers <- c("HBZ", "AHSP", "HBA1", "ALAS2", "CELF2", "FLI1")

p_markers <- FeaturePlot(combined_joined, features = markers, reduction = "umap",
                         pt.size = 0.3, ncol = 3)
p_markers

combined_joined$SingleR_labels_extended <- combined_joined$SingleR_labels_other

# overwrite for cluster 7
combined_joined$SingleR_labels_extended[Idents(combined_joined) == 7] <- "blood_progenitor_cell"

# plot using the new column
DimPlot(
  combined_joined,
  reduction = "umap",
  group.by  = "SingleR_labels_extended",
  label     = TRUE
)


p <- DimPlot(
  combined_joined,
  reduction = "umap",
  group.by = c("SingleR_labels_extended", "seurat_clusters"),
  label = TRUE
)
p


all.clusters.diff.genes <- FindAllMarkers(combined_joined, assay = "RNA", future.seed = FALSE)


# Top 10 markers per cluster by fold change (the genes that DEFINE each cluster).
all.clusters.diff.genes %>%
  group_by(cluster) %>%
  top_n(wt = avg_log2FC, n = 10)

# Top 10 by pct.1 = fraction of cells in the cluster expressing the gene
# (useful for finding clean, highly-detected markers).
all.clusters.diff.genes %>%
  group_by(cluster) %>%
  top_n(wt = pct.1, n = 10)

