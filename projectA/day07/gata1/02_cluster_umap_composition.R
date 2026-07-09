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
# -----------------------------------------------------------------------------
# HOW TO USE THIS TEMPLATE
#   Work through it top to bottom alongside the worksheet
#   (02_cluster_umap_composition.md). The boilerplate is already written for
#   you. Wherever you see a block like:
#
#       # ---- Step 1a: normalize ... ----
#       # Hint: NormalizeData(); args ...
#       # YOUR CODE HERE:
#       combined <-
#
#   finish the line yourself. If you get stuck, the completed answer key is in
#   scripts_finished/02_cluster_umap_composition.R -- try it on your own first.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

#library(future)
#plan("multisession", workers = 4)
# Also ensure BLAS sees the cores:
#Sys.setenv(OMP_NUM_THREADS = 4, OPENBLAS_NUM_THREADS = 4, MKL_NUM_THREADS = 4)
#options(future.globals.maxSize = 50 * 1024^3)  # 50 GB


library(Seurat)
library(dplyr)
library(ggplot2)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_qc.rds"))
DefaultAssay(combined) <- "RNA"

# ---- 1. Normalize -> variable genes -> scale --------------------------------
# (NormalizeData was already run in script 01, but re-running is harmless and
#  keeps this script runnable on its own.)

# ---- Step 1a: normalize (LogNormalize, scale factor 1e4) --------------------
# Cells differ in total counts for technical reasons; normalizing makes
# "expression" comparable across cells.
# Hint: NormalizeData(); args normalization.method = "LogNormalize", scale.factor = 1e4
# YOUR CODE HERE:
combined <-

# ---- Step 1b: highly variable genes (vst, 3000) -----------------------------
# Most genes are uninformative housekeepers; we model the 3000 most variable,
# where the biology lives.
# Hint: FindVariableFeatures(); args selection.method = "vst", nfeatures = 3000
# YOUR CODE HERE:
combined <-

head(VariableFeatures(combined), 20)

# ---- Step 1c: cell-cycle scoring --------------------------------------------
# Dividing cells can form their own clusters that have nothing to do with the
# biology we care about. Score each cell for S-phase and G2M-phase so we can
# see (and, if needed, regress out) that effect. Seurat ships curated gene
# lists; we hand them in below.
s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Hint: CellCycleScoring(); args s.features = s.genes, g2m.features = g2m.genes,
#       set.ident = TRUE. Adds S.Score, G2M.Score, and Phase to the metadata.
# YOUR CODE HERE:
combined <-

# ---- Step 1d: scale + regress out technical drivers -------------------------
# Center/scale each gene to mean 0, variance 1, and regress out sequencing
# depth and mitochondrial fraction so those technical drivers don't masquerade
# as biological structure in the PCA.
# Hint: ScaleData(); args features = VariableFeatures(combined),
#       vars.to.regress = c("nCount_RNA", "percent.mt")
# YOUR CODE HERE:
combined <-

# (Optional experiment: also regress out S.Score + G2M.Score to remove the
#  cell-cycle effect. Try it and compare your UMAPs.)
#combined <- ScaleData(
#  combined,
#  features = VariableFeatures(combined),
#  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
#  verbose = TRUE
#)

# ---- 2. PCA + how many PCs to keep ------------------------------------------
# ---- Step 2a: run PCA on the variable genes ---------------------------------
# Compress 3000 genes into a few dozen principal components.
# Hint: RunPCA(); args features = VariableFeatures(combined), npcs = 50, verbose = FALSE
# YOUR CODE HERE:
combined <-

print(combined[["pca"]], dims = 1:5, nfeatures = 10)

# Elbow plot: how much variance each PC explains. Where the curve flattens is
# roughly where added PCs become noise. (This diagnostic is written for you.)
ggsave(file.path(OUT_DIR, "gata1_elbow.png"),
       plot = ElbowPlot(combined, ndims = 50), width = 6, height = 4, dpi = 150)

plot = ElbowPlot(combined, ndims = 50)
plot

# ---- Step 2b: choose how many PCs to keep -----------------------------------
# The elbow flattens early here, so ~10 PCs capture most of the structure (20
# also works). The goal is that your story holds no matter which you pick.
# Hint: set n_pcs to a single number (e.g. 10).
# YOUR CODE HERE:
n_pcs <-

# ---- 3. UMAP + graph-based clustering ---------------------------------------
# ---- Step 3a: UMAP embedding ------------------------------------------------
# Project the n_pcs-dimensional PCA space into 2D for visualization.
# Hint: RunUMAP(); args dims = 1:n_pcs
# YOUR CODE HERE:
combined <-

# ---- Step 3b: build the neighbor graph --------------------------------------
# Hint: FindNeighbors(); args dims = 1:n_pcs
# YOUR CODE HERE:
combined <-

# ---- Step 3c: cluster the cells ---------------------------------------------
# Higher resolution -> more, smaller clusters. 0.4 is a moderate starting point.
# Hint: FindClusters(); args resolution = 0.4
# YOUR CODE HERE:
combined <-
#cluster numbers are random,
#so each time you run the FindClusters command above the groups will be about the same
#but which cluster number is which group will change



# Now lets look at the data!
#group.by uses the metadata we parsed in script 01.
#reduction can be set to pca or umap
colnames(combined@meta.data)

# Small helper (written for you): saves a plot to OUT_DIR at a fixed size.
save_dim <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
}

# ---- Step 3d: UMAP colored by cluster (WORKED EXAMPLE) ----------------------
# This first DimPlot is done for you. Notice the pattern:
#   DimPlot(object, reduction = "umap", group.by = <a metadata column>)
# then print the plot and save it. Copy this pattern for the ones below.
p <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p
save_dim(p, "gata1_umap_clusters.png")

# ---- Step 3e: UMAP colored by day -------------------------------------------
# Hint: same DimPlot pattern, group.by = "day"
# YOUR CODE HERE:
p <-
p
save_dim(p, "gata1_umap_day.png")

# ---- Step 3f: UMAP colored by construct -------------------------------------
# Hint: group.by = "construct"
# YOUR CODE HERE:
p <-
p
save_dim(p, "gata1_umap_construct.png")

# ---- Step 3g: UMAP colored by genotype --------------------------------------
# Hint: group.by = "genotype"
# YOUR CODE HERE:
p <-
p
save_dim(p, "gata1_umap_genotype.png")

# ---- Step 3h: cell-cycle score scatter (written for you) --------------------
# A quick look at S.Score vs G2M.Score to see how much cycling is going on.
p <- ggplot(combined@meta.data, aes(x = G2M.Score, y = S.Score)) + geom_hex()
p
save_dim(p, "gata1_umap_ccscatter.png")

# ---- Step 3i: UMAP colored by cell-cycle Phase ------------------------------
# Are the clusters driven by cell cycle? Color by Phase to check.
# Hint: group.by = "Phase"
# YOUR CODE HERE:
p <-
p
save_dim(p, "gata1_umap_ccPhase.png")

#this one takes a while because its not a factor
p <- DimPlot(combined, reduction = "umap", group.by = "ApopScore1")
p
save_dim(p, "gata1_umap_ApopScore1.png")

# ---- Step 3j: UMAP of day, split by construct_genotype ----------------------
# Facet the day-colored UMAP into one panel per construct_genotype so you can
# compare conditions side by side.
# Hint: DimPlot(..., group.by = "day", split.by = "construct_genotype",
#               label = TRUE, pt.size = 0.3)
# YOUR CODE HERE:
p <-
p
save_dim(p, "gata1_umap_day_by_construct_genotype.png", w = 12, h = 5)

# ---- 4. Composition: who lives in which cluster? ----------------------------
# Cross-tabulate cluster x sample. Then normalize two different ways, because
# the two normalizations answer two different questions.
# ---- Step 4a: cross-tabulate cluster x sample -------------------------------
# Hint: table(); rows = Idents(combined), columns = combined$sample
# YOUR CODE HERE:
ct <-

# now create a dataframe that divides by the proportion in each row
# (a) margin = 1 -> rows (clusters) sum to 1.
#     "Within this cluster, what fraction comes from each sample?"
# ---- Step 4b: normalize by cluster (rows sum to 1) --------------------------
# Hint: prop.table(); args ct, margin = 1
# YOUR CODE HERE:
frac_by_cluster <-

# ---- Step 4c: turn the table into a tidy data frame -------------------------
# Hint: as.data.frame(frac_by_cluster), then set colnames to
#       c("cluster", "sample", "fraction")
# YOUR CODE HERE:
frac_df <-
colnames(frac_df) <-

# Order the samples by similarity using hierarchical clustering, so visually
# similar columns sit next to each other instead of in arbitrary order.
# (The distance + hclust call is written for you.)
m <- as.matrix(frac_by_cluster)
hc_col    <- hclust(dist(t(m), method = "euclidean"), method = "complete")
col_order <- hc_col$labels[hc_col$order]

# ---- Step 4d: reorder the axes by the dendrogram ----------------------------
# Turn sample/cluster into factors whose level order follows the clustering.
# Hint: factor(frac_df$sample, levels = col_order) for the columns;
#       factor(frac_df$cluster, levels = sort(unique(frac_df$cluster))) for rows
# YOUR CODE HERE:
frac_df$sample  <-
frac_df$cluster <-

p_cluster <- ggplot(frac_df, aes(x = sample, y = cluster, fill = fraction)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Fraction of each cluster contributed by each sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_cluster
save_dim(p_cluster, "gata1_composition_by_cluster.png", w = 9, h = 6)

# now create a dataframe that divides by the proportion in each column
# (b) margin = 2 -> columns (samples) sum to 1.
#     "Within this sample, what fraction of cells fall in each cluster?"
#     This is the one that reveals shifts in cell-type proportions across
#     genotype / construct / day.
# ---- Step 4e: normalize by sample (columns sum to 1) ------------------------
# Hint: prop.table(); args ct, margin = 2. Then as.data.frame() and set
#       colnames to c("cluster", "sample", "fraction_sample").
# YOUR CODE HERE:
frac_by_sample <-
frac_df_sample <-
colnames(frac_df_sample) <-

# Order both axes by hierarchical clustering (written for you).
m2 <- as.matrix(frac_by_sample)
row_order <- { hc <- hclust(dist(m2),      method = "complete"); hc$labels[hc$order] }
col_order <- { hc <- hclust(dist(t(m2)),   method = "complete"); hc$labels[hc$order] }

# ---- Step 4f: reorder the axes by the dendrogram ----------------------------
# Hint: factor(frac_df_sample$cluster, levels = row_order) and
#       factor(frac_df_sample$sample,  levels = col_order)
# YOUR CODE HERE:
frac_df_sample$cluster <-
frac_df_sample$sample  <-

p_sample <- ggplot(frac_df_sample, aes(x = sample, y = cluster, fill = fraction_sample)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = "Cluster proportions within each sample (columns sum to 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_sample
save_dim(p_sample, "gata1_composition_by_sample.png", w = 9, h = 6)

# ---- 5. Save for the annotation script --------------------------------------
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_clustered.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_clustered.rds"))
