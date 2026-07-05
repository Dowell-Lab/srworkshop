# =============================================================================
# ADVANCED track - Script 06: Harmony integration + LISI/cLISI diagnostics
# -----------------------------------------------------------------------------
# WHAT THIS TEACHES
#   A different integration strategy from the FBM basic track (which used RPCA).
#   Harmony corrects batch effects in PCA space, and - crucially - this script
#   shows how to QUANTIFY whether integration worked using LISI:
#     * iLISI (integration LISI over the batch variable): higher = batches mixed
#     * cLISI (cell-type LISI over the label variable):   lower  = cell types kept distinct
#   We sweep a range of Harmony `theta` values and compare, so integration is a
#   measured decision, not a vibe.
#
# DATASET: fetal bone marrow (T21 + D21) pre-labeled objects. These come from
#   different conditions/donors, so there is a real batch effect worth removing
#   while preserving the published cell-type labels.
#
# STYLE: follows the lab's parameter-block ("io") pattern, with small helper
#   functions and a per-theta output loop, adapted to the course conventions
#   (source ../00_paths_and_setup.R, read-only inputs, write to OUT_DIR).
# =============================================================================

source("../00_paths_and_setup.R")

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(data.table)
library(RColorBrewer)
# LISI: install once with remotes::install_github("immunogenomics/lisi")
library(lisi)

# ---- Parameter block --------------------------------------------------------
# Everything you might tune lives here, mirroring the lab's io-list convention.
io <- list(
  out_dir              = OUT_DIR,
  sample_id            = "sample",        # per-cell batch/sample label we create below
  ident_col            = "cell.labels",   # published cell-type label on the objects
  integration_variable = "sample",        # variable Harmony mixes over (the "batch")
  contrast             = "condition",     # biological contrast (T21 vs D21)
  npcs                 = 20,
  resolution           = 0.4,
  thetas               = c(0, 1, 2, 3, 5),# Harmony diversity-penalty sweep
  harmony_seed         = 42L,
  marker_genes         = c("CD34", "GATA1", "HBB", "HBG1", "MPO", "CD79A",
                           "PTPRC", "GYPA", "PF4", "IRF8")
)
outdir_harmony <- file.path(io$out_dir, "advanced_harmony_integration")
dir.create(outdir_harmony, showWarnings = FALSE, recursive = TRUE)

# ---- Small helpers ----------------------------------------------------------
create_outdir <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

# UMAP + clustering off a chosen reduction (pca or harmony)
umap_cluster <- function(sobj, reduction, npcs, resolution) {
  sobj <- RunUMAP(sobj, reduction = reduction, dims = 1:npcs)
  sobj <- FindNeighbors(sobj, reduction = reduction, dims = 1:npcs)
  sobj <- FindClusters(sobj, resolution = resolution)
  sobj
}

# Compute per-cell LISI on the UMAP embedding for the requested labels
calculate_lisi <- function(sobj, cols) {
  scores <- compute_lisi(
    X              = Embeddings(sobj, "umap"),
    meta_data      = sobj@meta.data,
    label_colnames = cols
  )
  colnames(scores) <- paste0(colnames(scores), "_lisi")
  cbind(cellID = colnames(sobj), scores)
}

# ---- 1. Build ONE merged object from the two labeled FBM objects ------------
# load_one() (from 00_paths_and_setup.R) returns whatever single object each
# .RData file holds, so we're robust to internal naming.
t21 <- load_one(T21_LABELED)
d21 <- load_one(D21_LABELED)

# Stamp batch + biological condition so Harmony and LISI have something to read.
t21$sample <- "T21";  t21$condition <- "T21"
d21$sample <- "D21";  d21$condition <- "D21"

# Make sure the published label column is called ident_col on both.
# (The objects carry cell.labels / broad_extfig7A_cell.labels; adjust if needed.)
if (!io$ident_col %in% colnames(t21@meta.data)) {
  stop("Expected label column '", io$ident_col, "' not found on T21 object; ",
       "inspect colnames(t21@meta.data) and update io$ident_col.")
}

sobj_combo <- merge(t21, d21, add.cell.ids = c("T21", "D21"))
sobj_combo <- JoinLayers(sobj_combo)

# ---- 2. Standard preprocessing to get a PCA to hand to Harmony --------------
sobj_combo <- NormalizeData(sobj_combo, verbose = FALSE)
sobj_combo <- FindVariableFeatures(sobj_combo, nfeatures = 2000, verbose = FALSE)
sobj_combo <- ScaleData(sobj_combo, verbose = FALSE)
sobj_combo <- RunPCA(sobj_combo, npcs = io$npcs, verbose = FALSE)

# ---- 3. Loop: baseline (no integration) then each Harmony theta -------------
set.seed(io$harmony_seed)
subdirs  <- c("baseline", paste0("theta_", io$thetas))
lisi_all <- list()

for (outdir in subdirs) {
  cat("\n==== ", outdir, " ====\n")
  outdir_subdir <- file.path(outdir_harmony, outdir)
  create_outdir(outdir_subdir)

  if (outdir == "baseline") {
    # No correction: cluster/UMAP straight off PCA to see the batch effect.
    reduction_use <- "pca"
    sobj_combo <- umap_cluster(sobj_combo, "pca", io$npcs, io$resolution)
  } else {
    theta_val <- as.numeric(gsub("theta_", "", outdir))
    cat("Running Harmony with theta =", theta_val, "\n")
    sobj_combo <- RunHarmony(
      sobj_combo,
      group.by.vars    = io$integration_variable,
      theta            = theta_val,
      plot_convergence = TRUE,
      verbose          = TRUE
    )
    ggsave(file.path(outdir_subdir, "harmony_convergence.png"),
           width = 6, height = 5, dpi = 150)
    reduction_use <- "harmony"
    sobj_combo <- umap_cluster(sobj_combo, "harmony", io$npcs, io$resolution)
  }

  # UMAPs: colored by batch (should mix after integration) and by cell type
  # (should stay separated).
  ggsave(file.path(outdir_subdir, "umap_by_sample.png"),
         DimPlot(sobj_combo, group.by = io$sample_id) +
           ggtitle(paste0(outdir, " - by batch")),
         width = 6, height = 5, dpi = 150)
  ggsave(file.path(outdir_subdir, "umap_by_celltype.png"),
         DimPlot(sobj_combo, group.by = io$ident_col, label = TRUE, repel = TRUE) +
           NoLegend() + ggtitle(paste0(outdir, " - by cell type")),
         width = 7, height = 6, dpi = 150)
  ggsave(file.path(outdir_subdir, "umap_by_condition.png"),
         DimPlot(sobj_combo, group.by = io$contrast),
         width = 6, height = 5, dpi = 150)

  # ---- LISI: the quantitative integration check ----
  # iLISI over the batch variable  -> want HIGH (batches interspersed)
  # cLISI over the cell-type label -> want LOW  (types stay pure)
  scores <- as.data.table(calculate_lisi(sobj_combo,
                                          c(io$integration_variable, io$ident_col)))
  scores$stage <- outdir
  lisi_all[[outdir]] <- scores

  saveRDS(sobj_combo,
          file.path(outdir_subdir, paste0("sobj_", outdir, ".rds")))
}

# ---- 4. Summarize LISI across the theta sweep -------------------------------
# One tidy table + a summary plot so you can PICK a theta on evidence.
lisi_dt <- rbindlist(lisi_all, fill = TRUE)
ilisi_col <- paste0(io$integration_variable, "_lisi")
clisi_col <- paste0(io$ident_col, "_lisi")

lisi_summary <- lisi_dt[, .(
  median_iLISI = median(get(ilisi_col), na.rm = TRUE),
  median_cLISI = median(get(clisi_col), na.rm = TRUE)
), by = stage]
lisi_summary$stage <- factor(lisi_summary$stage, levels = subdirs)
print(lisi_summary)
fwrite(lisi_summary, file.path(outdir_harmony, "lisi_summary_by_theta.csv"))

# The sweet spot maximizes iLISI (mixing) while keeping cLISI near baseline.
p_sum <- ggplot(lisi_summary, aes(x = median_iLISI, y = median_cLISI, label = stage)) +
  geom_point(size = 3, color = "steelblue") +
  ggrepel::geom_text_repel() +
  labs(title = "Harmony theta sweep: batch mixing vs. cell-type purity",
       x = "median iLISI (batch mixing - higher is better)",
       y = "median cLISI (cell-type mixing - lower is better)") +
  theme_bw()
ggsave(file.path(outdir_harmony, "lisi_tradeoff.png"), p_sum, width = 7, height = 5, dpi = 150)

cat("\nFinished Harmony integration + LISI sweep.\n")
sessionInfo()
