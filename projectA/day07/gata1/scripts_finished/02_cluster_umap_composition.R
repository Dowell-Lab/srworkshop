# =============================================================================
# GATA1 track - Script 02: Normalize, PCA, cluster, UMAP, composition analysis
# -----------------------------------------------------------------------------
# Picks up the filtered object from script 01 and turns it into a clustered,
# UMAP-embedded dataset. Then it quantifies HOW each experimental condition
# is distributed across the clusters using composition heatmaps with
# hierarchical clustering.
#
# Inputs : OUT_DIR/gata1_combined_qc.rds  (from 01_load_qc_metadata.R)
# Outputs: a processed object + several plots in OUT_DIR.
# =============================================================================

source("../00_paths_and_setup.R")

library(Seurat)
library(dplyr)
library(ggplot2)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_qc.rds"))
DefaultAssay(combined) <- "RNA"

# ---- 1. Normalize -> variable genes -> scale --------------------------------
# (NormalizeData was already run in script 01, but re-running is harmless and
#  keeps this script runnable on its own.)
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)

# Highly variable genes carry most of the biological signal; we model 3000.
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                 nfeatures = 3000, verbose = FALSE)
head(VariableFeatures(combined), 20)

# Scale + regress out technical drivers so they don't dominate the PCA.
combined <- ScaleData(combined,
                      vars.to.regress = c("nCount_RNA", "percent.mt"),
                      verbose = FALSE)

# ---- 2. PCA + how many PCs to keep ------------------------------------------
combined <- RunPCA(combined, features = VariableFeatures(combined),
                   npcs = 50, verbose = FALSE)

print(combined[["pca"]], dims = 1:5, nfeatures = 10)
ggsave(file.path(OUT_DIR, "gata1_elbow.png"),
       plot = ElbowPlot(combined, ndims = 50), width = 6, height = 4, dpi = 150)

# The elbow flattens early here, so 10 PCs capture the structure.
n_pcs <- 10

# ---- 3. UMAP + graph-based clustering ---------------------------------------
combined <- RunUMAP(combined, dims = 1:n_pcs)
combined <- FindNeighbors(combined, dims = 1:n_pcs)
combined <- FindClusters(combined, resolution = 0.4)

# A few orientation plots. group.by uses the metadata we parsed in script 01.
save_dim <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
}
save_dim(DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE),
         "gata1_umap_clusters.png")
save_dim(DimPlot(combined, reduction = "umap", group.by = "day"),       "gata1_umap_day.png")
save_dim(DimPlot(combined, reduction = "umap", group.by = "construct"), "gata1_umap_construct.png")
save_dim(DimPlot(combined, reduction = "umap", group.by = "genotype"),  "gata1_umap_genotype.png")
save_dim(DimPlot(combined, reduction = "umap", group.by = "day", split.by = "construct_genotype",
                 label = TRUE, pt.size = 0.3),
         "gata1_umap_day_by_construct_genotype.png", w = 12, h = 5)

# ---- 4. Composition: who lives in which cluster? ----------------------------
# Cross-tabulate cluster x sample. Then normalize two different ways, because
# the two normalizations answer two different questions.
ct <- table(Idents(combined), combined$sample)

# (a) margin = 1 -> rows (clusters) sum to 1.
#     "Within this cluster, what fraction comes from each sample?"
frac_by_cluster <- prop.table(ct, margin = 1)
frac_df <- as.data.frame(frac_by_cluster)
colnames(frac_df) <- c("cluster", "sample", "fraction")

# Order the samples by similarity using hierarchical clustering, so visually
# similar columns sit next to each other instead of in arbitrary order.
m <- as.matrix(frac_by_cluster)
hc_col    <- hclust(dist(t(m), method = "euclidean"), method = "complete")
col_order <- hc_col$labels[hc_col$order]

frac_df$sample  <- factor(frac_df$sample,  levels = col_order)
frac_df$cluster <- factor(frac_df$cluster, levels = sort(unique(frac_df$cluster)))

p_cluster <- ggplot(frac_df, aes(x = sample, y = cluster, fill = fraction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Fraction of each cluster contributed by each sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_dim(p_cluster, "gata1_composition_by_cluster.png", w = 9, h = 6)

# (b) margin = 2 -> columns (samples) sum to 1.
#     "Within this sample, what fraction of cells fall in each cluster?"
#     This is the one that reveals shifts in cell-type proportions across
#     genotype / construct / day.
frac_by_sample <- prop.table(ct, margin = 2)
frac_df_sample <- as.data.frame(frac_by_sample)
colnames(frac_df_sample) <- c("cluster", "sample", "fraction_sample")

m2 <- as.matrix(frac_by_sample)
row_order <- { hc <- hclust(dist(m2),      method = "complete"); hc$labels[hc$order] }
col_order <- { hc <- hclust(dist(t(m2)),   method = "complete"); hc$labels[hc$order] }

frac_df_sample$cluster <- factor(frac_df_sample$cluster, levels = row_order)
frac_df_sample$sample  <- factor(frac_df_sample$sample,  levels = col_order)

p_sample <- ggplot(frac_df_sample, aes(x = sample, y = cluster, fill = fraction_sample)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = "Cluster proportions within each sample (columns sum to 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_dim(p_sample, "gata1_composition_by_sample.png", w = 9, h = 6)

# ---- 5. Save for the annotation script --------------------------------------
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_clustered.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_clustered.rds"))
