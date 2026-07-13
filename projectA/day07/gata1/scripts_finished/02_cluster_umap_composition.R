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
#
# ANSWER KEY for the fill-in template ../02_cluster_umap_composition.R. Try the
# template (with the worksheet 02_cluster_umap_composition.md) before reading
# this.
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

# If you just finished script 01 in this same R session, `combined` is already
# in your environment -- you do NOT need to reload it, so comment out the next
# line. Only run it on a fresh session (or if you cleared your workspace).
combined <- readRDS(file.path(OUT_DIR, "gata1_combined_qc.rds"))

#If you didn't make this yet, use mine!
#combined <- readRDS(file.path(COOKING, "gata1_combined_qc.rds"))


DefaultAssay(combined) <- "RNA"

# Join the per-sample layers (Seurat v5) so the whole workflow below runs on one
# combined layer. Harmless no-op if they are already joined.
combined <- JoinLayers(combined)

# ---- 1. Normalize -> variable genes -> scale --------------------------------
# (NormalizeData was already run in script 01, but re-running is harmless and
#  keeps this script runnable on its own.)
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = TRUE)

# Highly variable genes carry most of the biological signal; we model 3000.
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                 nfeatures = 3000, verbose = TRUE)
head(VariableFeatures(combined), 20)

s.genes  <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

combined <- CellCycleScoring(
  combined,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# Scale + regress out technical drivers so they don't dominate the PCA.
# HEADS UP: with the regression this is the SLOW step -- it can take a while
# (up to ~20 minutes). verbose = TRUE lets you watch its progress.
combined <- ScaleData(
  combined,
  features = VariableFeatures(combined),
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = TRUE
)

#combined <- ScaleData(
#  combined,
#  features = VariableFeatures(combined),
#  vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
#  verbose = TRUE
#)

# ---- 2. PCA + how many PCs to keep ------------------------------------------
# HEADS UP: PCA also takes a while on this many cells. verbose = TRUE prints the
# top genes for each PC as it goes, so you can see it working.
combined <- RunPCA(combined, features = VariableFeatures(combined),
                   npcs = 50, verbose = TRUE)

print(combined[["pca"]], dims = 1:5, nfeatures = 10)

# We just VIEW the elbow below; uncomment the ggsave to write gata1_elbow.png.
# ggsave(file.path(OUT_DIR, "gata1_elbow.png"),
#        plot = ElbowPlot(combined, ndims = 50), width = 6, height = 4, dpi = 150)

p_elbow <- ElbowPlot(combined, ndims = 50)
p_elbow

# The elbow flattens early here, so 10 PCs capture most of the structure, you can also use 20
# the goal is that your story holds no matter what n_pcs you choose
n_pcs <- 10

# ---- 3. UMAP + graph-based clustering ---------------------------------------
combined <- RunUMAP(combined, dims = 1:n_pcs)
combined <- FindNeighbors(combined, dims = 1:n_pcs)
combined <- FindClusters(combined, resolution = 0.4)
#cluster numbers are random, 
#so each time you run the FindClusters command above the groups will be about the same
#but which cluster number is which group will change



# Now lets look at the data!
#group.by uses the metadata we parsed in script 01.
#reduction can be set to pca or umap
colnames(combined@meta.data)

# For the workshop we just VIEW each plot (the bare `p` line above every
# save_dim() call prints it), so the ggsave is disabled -- save_dim() is a no-op
# here. Uncomment the ggsave line to write the PNGs to OUT_DIR instead.
save_dim <- function(p, file, w = 7, h = 5) {
  # ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
  invisible(p)
}

p<- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p
save_dim(p,"gata1_umap_clusters.png")

p<-DimPlot(combined, reduction = "umap", group.by = "day")
p
save_dim(p,"gata1_umap_day.png")

p<- DimPlot(combined, reduction = "umap", group.by = "construct")
p
save_dim(p,"gata1_umap_construct.png")


p<- DimPlot(combined, reduction = "umap", group.by = "genotype")
p
save_dim(p,"gata1_umap_genotype.png")


p <- ggplot(combined@meta.data, aes(x=G2M.Score, y=S.Score))+geom_hex()
p

save_dim(p,"gata1_umap_ccscatter.png")



p<- DimPlot(combined, reduction = "umap", group.by = "Phase")
p
save_dim(p,"gata1_umap_ccPhase.png")

#this one takes a while because its not a factor
#p<- DimPlot(combined, reduction = "umap", group.by = "ApopScore1")
#p
#save_dim(p,"gata1_umap_ApopScore1.png")

p<- DimPlot(combined, reduction = "umap", group.by = "day", split.by = "construct_genotype",
            label = TRUE, pt.size = 0.3)
p
save_dim(p,"gata1_umap_day_by_construct_genotype.png", w = 12, h = 5)

# ---- 4. Composition: who lives in which cluster? ----------------------------
# Cross-tabulate cluster x sample. Then normalize two different ways, because
# the two normalizations answer two different questions.
ct <- table(Idents(combined), combined$sample)

# now create a dataframe that divides by the proportion in each row
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
p_cluster
save_dim(p_cluster, "gata1_composition_by_cluster.png", w = 9, h = 6)

# now create a dataframe that divides by the proportion in each column
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
p_sample
save_dim(p_sample, "gata1_composition_by_sample.png", w = 9, h = 6)

# ---- 5. Save for the annotation script --------------------------------------
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_clustered.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_clustered.rds"))
