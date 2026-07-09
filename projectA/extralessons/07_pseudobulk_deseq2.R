# =============================================================================
# ADVANCED track - Script 07: Pseudobulk differential expression with DESeq2
# -----------------------------------------------------------------------------
# WHAT THIS TEACHES
#   The statistically correct way to compare CONDITIONS in single-cell data:
#   collapse cells to one count profile per (cell type x sample) - "pseudobulk" -
#   then run a replicate-aware negative-binomial test with DESeq2.
#
#   Why not just use FindMarkers across conditions? Because thousands of cells
#   from ONE sample are not thousands of independent replicates. Treating them
#   as such inflates significance massively. Pseudobulk restores the real unit
#   of replication (the sample) and controls false positives.
#   (Squair et al., Nat Commun 2021: https://www.nature.com/articles/s41467-021-25960-2)
#
# DATASET: GATA1 (GSE271399). Its clean factorial design
#   genotype x construct x day gives real biological replicates to aggregate.
#   Example contrast here: GATA1s vs wtGATA1 (the "construct" effect).
#
# INPUT : OUT_DIR/gata1_combined_annotated.rds  (from gata1/03)
# STYLE : follows the lab's io-block + edgeR::filterByExpr + rlog QC + ashr
#         shrinkage + fdrtool pattern, adapted to course conventions.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(data.table)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
# ashr (lfcShrink type="ashr") and fdrtool are optional but recommended:
#   install.packages(c("ashr", "fdrtool"))

# ---- Parameter block --------------------------------------------------------
io <- list(
  sobj_rds     = file.path(OUT_DIR, "gata1_combined_annotated.rds"),
  out_dir      = OUT_DIR,
  sample_id    = "sample",                # pseudobulk replicate unit
  ident_col    = "SingleR_labels_other",  # cell type to test WITHIN
  contrast     = "construct",             # factor being compared
  de_model     = "construct",             # DESeq2 design (add covariates e.g. "day + construct")
  assay        = "RNA",
  min_count    = 10,                      # edgeR filterByExpr threshold
  test_level   = "GATA1s",                # numerator of the contrast
  ref_level    = "wtGATA1"                # denominator (reference)
)
outdir_pseudo <- file.path(io$out_dir,
                           paste0("advanced_pseudobulk_", io$ident_col, "_", io$de_model))
dir.create(outdir_pseudo, showWarnings = FALSE, recursive = TRUE)

save_pheatmap_pdf <- function(x, filename, w = 7, h = 7) {
  pdf(filename, width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(x$gtable); dev.off()
}

# ---- 1. Load and prepare ----------------------------------------------------
sobj <- readRDS(io$sobj_rds)
sobj <- JoinLayers(sobj)
Idents(sobj) <- io$ident_col
DefaultAssay(sobj) <- io$assay

# Drop the "other" bucket so we test real cell types only.
sobj <- subset(sobj, idents = setdiff(unique(Idents(sobj)), c("other", NA)))

# ---- 2. Pseudobulk: sum counts per (cell type x sample x contrast) ----------
# AggregateExpression returns SUMMED raw counts - exactly what DESeq2 wants.
group.factors <- unique(c(io$ident_col, io$sample_id,
                          trimws(strsplit(io$de_model, "+", fixed = TRUE)[[1]])))

bulk <- AggregateExpression(
  sobj,
  assays        = io$assay,
  return.seurat = TRUE,
  group.by      = group.factors
)
# DESeq2 needs integers.
bulk[[io$assay]]$counts <- round(bulk[[io$assay]]$counts)

md <- as.data.table(bulk@meta.data, keep.rownames = "pseudobulk_id")
print(md[, .N, by = get(io$ident_col)])

# ---- 3. Loop over cell types; run DESeq2 within each -------------------------
res_all <- list()
for (myclust in unique(md[[io$ident_col]])) {
  cat("\n==== Cell type:", myclust, "====\n")
  outdir_base <- file.path(outdir_pseudo, gsub("[^A-Za-z0-9]+", "_", myclust))
  dir.create(outdir_base, showWarnings = FALSE, recursive = TRUE)

  md.sub  <- md[get(io$ident_col) == myclust, ]
  # Need at least 2 replicates per side to fit the model.
  tab <- table(md.sub[[io$contrast]])
  if (length(tab) < 2 || any(tab < 2)) {
    cat("  Skipping - fewer than 2 replicates per group.\n"); next
  }
  counts <- bulk[[io$assay]]$counts[, md.sub$pseudobulk_id, drop = FALSE]

  # ---- 3a. Gene prefilter with edgeR (group-aware) ----
  grp  <- factor(md.sub[[io$contrast]])
  dge  <- DGEList(counts = counts, group = grp)
  keep <- filterByExpr(dge, group = grp, min.count = io$min_count)
  cat("  Keeping", sum(keep), "/", length(keep), "genes.\n")
  counts <- counts[keep, , drop = FALSE]

  # ---- 3b. Build DESeq2 object; set reference level explicitly ----
  md.sub[[io$contrast]] <- relevel(factor(md.sub[[io$contrast]]), ref = io$ref_level)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = as.data.frame(md.sub),
    design    = as.formula(paste0("~", io$de_model))
  )
  # poscounts size factors are robust to the many zeros in pseudobulk data.
  dds <- tryCatch(DESeq(dds),
                  error = function(e) {
                    message("  Default DESeq failed; retrying sfType='poscounts'.")
                    DESeq(dds, sfType = "poscounts")
                  })
  saveRDS(dds, file.path(outdir_base, "deseq2_dds.rds"))

  # dispersion diagnostic
  pdf(file.path(outdir_base, "dispersion.pdf"), width = 5, height = 4)
  plotDispEsts(dds); dev.off()

  # ---- 3c. Unsupervised QC on rlog counts (do samples group sensibly?) ----
  rld <- rlog(dds, blind = TRUE)
  pca_df <- plotPCA(rld, intgroup = io$contrast, returnData = TRUE)
  pv <- round(100 * attr(pca_df, "percentVar"))
  ggsave(file.path(outdir_base, "rlog_PCA.pdf"),
         ggplot(pca_df, aes(PC1, PC2, color = group)) + geom_point(size = 4) +
           labs(title = myclust, x = paste0("PC1 ", pv[1], "%"),
                y = paste0("PC2 ", pv[2], "%")) + theme_classic(),
         width = 5, height = 4)
  d <- dist(t(assay(rld)))
  ph <- pheatmap(as.matrix(d), main = paste0("Sample distances: ", myclust),
                 col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                 clustering_method = "ward.D")
  save_pheatmap_pdf(ph, file.path(outdir_base, "sample_distances.pdf"))

  # ---- 3d. Extract the contrast (GATA1s vs wtGATA1) ----
  res <- results(dds, contrast = c(io$contrast, io$test_level, io$ref_level))
  res.dt <- as.data.table(res, keep.rownames = "Gene")

  # p-value histogram - a healthy test has a flat null + spike near 0
  ggsave(file.path(outdir_base, "pvalue_histogram.pdf"),
         ggplot(res.dt, aes(pvalue)) +
           geom_histogram(bins = 30, fill = "skyblue", color = "black") +
           ggtitle(paste0(myclust, ": ", io$test_level, " vs ", io$ref_level)) +
           theme_bw(),
         width = 4, height = 3)

  # optional local-FDR recalibration (fdrtool)
  if (requireNamespace("fdrtool", quietly = TRUE)) {
    ok <- !is.na(res.dt$pvalue)
    fdr <- fdrtool::fdrtool(res.dt$pvalue[ok], statistic = "pvalue", plot = FALSE)
    res.dt[ok, c("qvalue", "lfdr") := list(fdr$qval, fdr$lfdr)]
  }

  # ---- 3e. LFC shrinkage (ashr) for stable effect sizes ----
  if (requireNamespace("ashr", quietly = TRUE)) {
    res_lfc <- lfcShrink(dds, contrast = c(io$contrast, io$test_level, io$ref_level),
                         type = "ashr")
    res.dt$log2FoldChange_ashr <- res_lfc$log2FoldChange
  }

  # MA plot
  toplot <- as.data.frame(res.dt)
  toplot$DEG <- ifelse(!is.na(toplot$padj) & toplot$padj < 0.05, "padj<0.05", "NS")
  ggsave(file.path(outdir_base, "MA_plot.pdf"),
         ggplot(toplot, aes(log2(baseMean), log2FoldChange, color = DEG)) +
           geom_hline(yintercept = 0, lty = 2, color = "grey") +
           geom_point(alpha = 0.6) +
           scale_color_manual(values = c("NS" = "grey60", "padj<0.05" = "red3")) +
           theme_classic(), width = 5, height = 4)

  fwrite(res.dt, file.path(outdir_base, "deseq2_results.csv"))
  res_all[[myclust]] <- cbind(cell_type = myclust, res.dt)
}

# ---- 4. Merge results across cell types -------------------------------------
if (length(res_all) > 0) {
  merged <- rbindlist(res_all, fill = TRUE)
  fwrite(merged, file.path(outdir_pseudo, "deseq2_results_all_celltypes.csv"))
  cat("\nWrote merged results:", nrow(merged), "rows.\n")
}
cat("Finished pseudobulk DESeq2.\n")
sessionInfo()
