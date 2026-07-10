# =============================================================================
# GATA1 track - Script 01: Load, QC, metadata-from-names, and module scores
# -----------------------------------------------------------------------------
# Dataset: GSE271399 (Takasaki et al., Stem Cell Reports 2025)
#   Isogenic iPSC-derived hematopoietic progenitors.
#   Design = genotype (Euploid / T21) x construct (wtGATA1 / GATA1s)
#            x differentiation day (D7 / D9 / D11)  ->  12 conditions.
#
# This is the PRIMARY worked example for the course. You run it end to end.
#
# What this script does:
#   1. Reads each sample's 10x-style matrix from the READ-ONLY shared folder.
#   2. Builds one Seurat object per sample, then merges them.
#   3. Computes mitochondrial content and QC violin plots.
#   4. Filters cells on nFeature_RNA and percent.mt.
#   5. Derives clean metadata (genotype / construct / day) FROM the sample
#      names with code, instead of typing a spreadsheet by hand.
#   6. Adds stress and apoptosis module scores to flag a low-quality sample.
#
# Nothing is downloaded; the GSE271399 files already live on the cluster.
# All outputs are written to OUT_DIR (your own writable folder).
#
# ANSWER KEY for the fill-in template ../01_load_qc_metadata.R. Try the template
# (with the worksheet 01_load_qc_metadata.md) before reading this.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- 1. Read every sample from the shared, read-only folder ----------------
# Files are named GSE271399_<sample>_{matrix.mtx,barcodes.tsv,features.tsv}.gz
# ReadMtx wants the three pieces explicitly. We loop over GATA1_SAMPLES
# (defined once in 00_paths_and_setup.R) so the names live in a single place.

#this is how you read in one file
mtx_file  <- file.path(GATA1_DIR, "GSE271399_EuploidGATA1sD7_matrix.mtx.gz")
bc_file   <- file.path(GATA1_DIR, "GSE271399_EuploidGATA1sD7_barcodes.tsv.gz")
feat_file <- file.path(GATA1_DIR, "GSE271399_EuploidGATA1sD7_features.tsv.gz")

mat <- ReadMtx(
  mtx            = mtx_file,
  cells          = bc_file,
  features       = feat_file,
  cell.column    = 1,   # barcodes are in column 1
  feature.column = 2,   # gene SYMBOL is column 2 in a 10x features.tsv
  mtx.transpose  = FALSE
)

obj <- CreateSeuratObject(counts = mat, project = "EuploidGATA1sD7")

#But that would be really slow?
#this is how you read in a lot of files with the same structure to the names

GATA1_SAMPLES

seurat_list <- lapply(GATA1_SAMPLES, function(smpl) {
  mtx_file  <- file.path(GATA1_DIR, paste0("GSE271399_", smpl, "_matrix.mtx.gz"))
  bc_file   <- file.path(GATA1_DIR, paste0("GSE271399_", smpl, "_barcodes.tsv.gz"))
  feat_file <- file.path(GATA1_DIR, paste0("GSE271399_", smpl, "_features.tsv.gz"))

  mat <- ReadMtx(
    mtx            = mtx_file,
    cells          = bc_file,
    features       = feat_file,
    cell.column    = 1,   # barcodes are in column 1
    feature.column = 2,   # gene SYMBOL is column 2 in a 10x features.tsv
    mtx.transpose  = FALSE
  )

  obj <- CreateSeuratObject(counts = mat, project = smpl)
  obj$sample <- smpl      # stamp the sample name onto every cell
  obj
})

names(seurat_list) <- GATA1_SAMPLES

seurat_list

# A merged object BEFORE filtering. We keep it so we can compare cell counts
# before vs. after QC. "ori" = original.
combinedori <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = GATA1_SAMPLES
)

# ---- 2. Mitochondrial percentage --------------------------------------------
# Human mito genes start with "MT-". A high fraction usually means a dying or
# lysed cell whose cytoplasmic RNA leaked out, leaving mostly mito transcripts.
mito_pattern <- "^MT-"

seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
  obj
})

# ---- 3. QC violin plots (one PNG per sample) --------------------------------
qc_dir <- file.path(OUT_DIR, "gata1_qc_violin_plots")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

#before you run this section you should make sure you can see plots
for (nm in names(seurat_list)) {
  p <- VlnPlot(
    seurat_list[[nm]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3, assay = "RNA", layer = "counts"
  ) + patchwork::plot_annotation(title = nm)

  ggsave(file.path(qc_dir, paste0("QC_violin_", nm, ".png")),
         plot = p, width = 8, height = 4, dpi = 150)
  print(p)
}

#something is strange about EuploidwtGATA1D7
#something is strange about T21GATA1sD11
#What is it?

#https://www.biostars.org/p/407036/

#Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. 
#High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet). 
#In combination with %mitochondrial reads, removing outliers from these groups removes most doublets/dead cells/empty droplets, hence why filtering is a common pre-processing step.


# ---- 4. Filter cells --------------------------------------------------------

#these are not numbers you can just use!!!
#You have to look at the QC_violin plots to decide.

# Keep cells with:
#   nFeature_RNA > 200   (drop empty droplets / debris) 
#   nFeature_RNA < 6000  (drop likely doublets)
#   percent.mt   < 15    (drop dying cells) 


qc_filtered <- lapply(seurat_list, function(obj) {
  subset(obj,
         subset = nFeature_RNA > 200 &
                  nFeature_RNA < 6000 &
                  percent.mt   < 15)
})
names(qc_filtered) <- GATA1_SAMPLES

combined <- merge(
  x            = qc_filtered[[1]],
  y            = qc_filtered[-1],
  add.cell.ids = names(qc_filtered)
)

# Seurat v5 keeps each sample's counts in a separate layer after merge. With a
# plain merge (no integration) we want them joined so NormalizeData /
# AddModuleScore below operate on one combined layer.
combined <- JoinLayers(combined)

# How many cells did each sample lose to QC?
before <- as.data.frame(table(combinedori$sample)); colnames(before) <- c("sample", "n_before")
after  <- as.data.frame(table(combined$sample));    colnames(after)  <- c("sample", "n_after")
qc_summary <- merge(before, after, by = "sample")
qc_summary$kept_pct <- round(100 * qc_summary$n_after / qc_summary$n_before, 1)
print(qc_summary)
write.csv(qc_summary, file.path(OUT_DIR, "gata1_qc_summary.csv"), row.names = FALSE)

# ---- 5. Metadata FROM the sample names (the clean way) ----------------------
# The sample name already encodes the full design, e.g. "T21wtGATA1D11".
# Instead of maintaining a hand-typed metadata sheet (which is how the messy
# fetal-bone-marrow dataset went wrong), we PARSE the design out of the name.
# One source of truth = the file name = impossible to mislabel.

combined$genotype <- ifelse(grepl("^Euploid", combined$sample), "Euploid", "T21")

combined$construct <- dplyr::case_when(
  grepl("wtGATA1", combined$sample) ~ "wtGATA1",
  grepl("GATA1s",  combined$sample) ~ "GATA1s",
  TRUE                              ~ NA_character_
)

combined$day <- dplyr::case_when(
  grepl("D7$",  combined$sample) ~ "D7",
  grepl("D9$",  combined$sample) ~ "D9",
  grepl("D11$", combined$sample) ~ "D11",
  TRUE                           ~ NA_character_
)

# Convenient combined factors for split.by / group.by later.
combined$day_construct          <- paste0(combined$day, "_", combined$construct)
combined$day_construct_genotype <- paste0(combined$day, "_", combined$construct, "_", combined$genotype)
combined$construct_genotype     <- paste0(combined$construct, "_", combined$genotype)

# Sanity check: every cell got a genotype/construct/day (no NAs).
stopifnot(!any(is.na(combined$genotype)))
stopifnot(!any(is.na(combined$construct)))
stopifnot(!any(is.na(combined$day)))
message("Metadata parsed cleanly for all ", ncol(combined), " cells.")

# ---- 6. Stress & apoptosis module scores ------------------------------------
# A module score = average expression of a gene set, minus a matched random
# background. It turns "many genes" into one interpretable number per cell.
# Here we use it to confirm that one sample (T21GATA1sD11) is genuinely
# stressed/dying rather than biologically interesting.
stress_genes    <- c("FOS", "JUN", "JUNB", "HSPA1A", "HSP90AA1", "ATF3", "EGR1")
apoptosis_genes <- c("BAX", "BAK1", "CASP3", "CASP8", "FAS", "BBC3")

# NormalizeData first so AddModuleScore reads sensible "data" values.
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)

combined <- AddModuleScore(combined, features = list(stress_genes),    name = "StressScore")
combined <- AddModuleScore(combined, features = list(apoptosis_genes), name = "ApopScore")

# Compare the suspect sample to everyone else.
p_stress <- VlnPlot(combined, features = c("StressScore1", "ApopScore1", "percent.mt"),
                    group.by = "sample", pt.size = 0) +
  patchwork::plot_annotation(title = "Stress / apoptosis / mito by sample")
ggsave(file.path(OUT_DIR, "gata1_module_scores_by_sample.png"),
       plot = p_stress, width = 12, height = 5, dpi = 150)



# ---- 7. Save the filtered, annotated object for the next script -------------
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_qc.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_qc.rds"))

