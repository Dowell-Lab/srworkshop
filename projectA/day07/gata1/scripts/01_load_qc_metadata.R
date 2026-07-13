# =============================================================================
# GATA1 track - Script 01: Load, QC, metadata-from-names, and module scores
# -----------------------------------------------------------------------------
# Dataset: GSE271399 (Takasaki et al., Stem Cell Reports 2025)
#   Isogenic iPSC-derived hematopoietic progenitors.
#   Design = genotype (Euploid / T21) x construct (wtGATA1 / GATA1s)
#            x differentiation day (D7 / D9 / D11)  ->  12 conditions.
#
# This is the PRIMARY worked example for the course. You run it end to end:
# it reads each sample's 10x matrix, merges them, computes QC, filters cells,
# parses clean metadata straight out of the sample names, and adds stress /
# apoptosis module scores. The saved object feeds script 02.
#
# Inputs : read-only 10x files under GATA1_DIR (nothing is downloaded).
# Outputs: OUT_DIR/gata1_combined_qc.rds  (+ QC violins and a summary CSV).
# -----------------------------------------------------------------------------
# HOW TO USE THIS TEMPLATE
#   Work through it top to bottom alongside the worksheet
#   (01_load_qc_metadata.md). The boilerplate is already written for you.
#   Wherever you see a block like:
#
#       # ---- Step 2a: mitochondrial percentage ... ----
#       # Hint: PercentageFeatureSet(); args ...
#       # YOUR CODE HERE:
#       obj[["percent.mt"]] <-
#
#   finish the line yourself. If you get stuck, the completed answer key is in
#   scripts_finished/01_load_qc_metadata.R -- try it on your own first.
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

# ---- Step 1a: read ONE sample (WORKED EXAMPLE) ------------------------------
# This first read is done for you. Notice the pattern: point ReadMtx at the
# three gzipped 10x files, then wrap the matrix in a Seurat object. The loop
# below repeats exactly this for every sample.
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

# ---- Read a lot of files with the same name structure -----------------------
# We loop over GATA1_SAMPLES, building the three paths for each sample. The
# ReadMtx call inside is the same one you just saw, written for you.
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

  # ---- Step 1b: build one Seurat object for this sample ---------------------
  # Wrap the counts matrix in a Seurat object, tagging it with the sample name.
  # Hint: CreateSeuratObject(); args counts = mat, project = smpl
  # YOUR CODE HERE:
  obj <-

  # ---- Step 1c: stamp the sample name onto every cell -----------------------
  # This single line is what makes every downstream group.by / split.by and the
  # metadata parsing in Step 5 possible. Without it, the merged object forgets
  # which cell came from which sample.
  # Hint: assign smpl into the object's "sample" metadata column: obj$sample
  # YOUR CODE HERE:
  obj$sample <-

  obj
})

names(seurat_list) <- GATA1_SAMPLES

seurat_list

# ---- Step 1d: merge all samples (before filtering) --------------------------
# A merged object BEFORE filtering. We keep it so we can compare cell counts
# before vs. after QC. "ori" = original. This is a plain merge with NO batch
# correction -- fine here because the lines are isogenic and processed together.
# Hint: merge(); args x = seurat_list[[1]], y = seurat_list[-1],
#       add.cell.ids = GATA1_SAMPLES
# YOUR CODE HERE:
combinedori <-

# ---- 2. Mitochondrial percentage --------------------------------------------
# Human mito genes start with "MT-". A high fraction usually means a dying or
# lysed cell whose cytoplasmic RNA leaked out, leaving mostly mito transcripts.
mito_pattern <- "^MT-"

# ---- Step 2a: add percent.mt to every sample --------------------------------
# Compute, per cell, the fraction of reads coming from mitochondrial genes and
# store it in the "percent.mt" metadata column.
# Hint: PercentageFeatureSet(); args obj, pattern = mito_pattern
# YOUR CODE HERE:
seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <-
  obj
})

# ---- 3. QC violin plots (one PNG per sample) --------------------------------
# These plots are written for you. LOOK at them per sample before you pick
# thresholds in Step 4 -- the spread is different for each sample.
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

#Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet.
#High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or multiplet).
#In combination with %mitochondrial reads, removing outliers from these groups
#removes most doublets/dead cells/empty droplets, hence why filtering is a common
#pre-processing step.

# ---- 4. Filter cells --------------------------------------------------------
# ---- Step 4a: subset each sample on QC thresholds ---------------------------
# these are NOT numbers you can just use!! Look at the QC_violin plots first.
# Reasonable STARTING values for these cultured progenitors:
#   nFeature_RNA > 200   (drop empty droplets / debris)
#   nFeature_RNA < 6000  (drop likely doublets)
#   percent.mt   < 15    (drop dying cells)
# Hint: subset(); args obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
#       percent.mt < 15
# YOUR CODE HERE:
qc_filtered <- lapply(seurat_list, function(obj) {

})
names(qc_filtered) <- GATA1_SAMPLES

# Re-merge the filtered objects (written for you; same pattern as Step 1d).
combined <- merge(
  x            = qc_filtered[[1]],
  y            = qc_filtered[-1],
  add.cell.ids = names(qc_filtered)
)

# Seurat v5 keeps each sample's counts in a separate layer after merge. With a
# plain merge (no integration) we want them joined so NormalizeData /
# AddModuleScore below operate on one combined layer. (written for you)
combined <- JoinLayers(combined)

# How many cells did each sample lose to QC? (written for you)
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

# ---- Step 5a: genotype ------------------------------------------------------
# Names starting with "Euploid" are Euploid; everything else is T21.
# Hint: ifelse(); args grepl("^Euploid", combined$sample), "Euploid", "T21"
# YOUR CODE HERE:
combined$genotype <-

# ---- Step 5b: construct -----------------------------------------------------
# The construct is either "wtGATA1" or "GATA1s", spelled out in the name.
# Hint: dplyr::case_when(); grepl("wtGATA1", ...) ~ "wtGATA1",
#       grepl("GATA1s", ...) ~ "GATA1s", TRUE ~ NA_character_
# YOUR CODE HERE:
combined$construct <-

# ---- Step 5c: differentiation day -------------------------------------------
# The day is the suffix: D7, D9, or D11. Anchor with "$" so "D7" doesn't also
# match inside another token.
# Hint: dplyr::case_when(); grepl("D7$", ...) ~ "D7", grepl("D9$", ...) ~ "D9",
#       grepl("D11$", ...) ~ "D11", TRUE ~ NA_character_
# YOUR CODE HERE:
combined$day <-

# Convenient combined factors for split.by / group.by later. (written for you)
combined$day_construct          <- paste0(combined$day, "_", combined$construct)
combined$day_construct_genotype <- paste0(combined$day, "_", combined$construct, "_", combined$genotype)
combined$construct_genotype     <- paste0(combined$construct, "_", combined$genotype)

# Sanity check: every cell got a genotype/construct/day (no NAs). These fail
# LOUDLY if a typo in the sample list left any cell unlabeled. (written for you)
stopifnot(!any(is.na(combined$genotype)))
stopifnot(!any(is.na(combined$construct)))
stopifnot(!any(is.na(combined$day)))
message("Metadata parsed cleanly for all ", ncol(combined), " cells.")

# ---- 6. Stress & apoptosis module scores ------------------------------------
# A module score = average expression of a gene set, minus a matched random
# background. It turns "many genes" into one interpretable number per cell.
# Here we use it to confirm that a suspect sample is genuinely stressed/dying
# rather than biologically interesting.
stress_genes    <- c("FOS", "JUN", "JUNB", "HSPA1A", "HSP90AA1", "ATF3", "EGR1")
apoptosis_genes <- c("BAX", "BAK1", "CASP3", "CASP8", "FAS", "BBC3")

# NormalizeData first so AddModuleScore reads sensible "data" values. (written for you)
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)

# ---- Step 6a: stress module score -------------------------------------------
# Score the immediate-early / heat-shock set that spikes during dissociation.
# Hint: AddModuleScore(); args combined, features = list(stress_genes),
#       name = "StressScore"  (adds a column "StressScore1")
# YOUR CODE HERE:
combined <-

# ---- Step 6b: apoptosis module score ----------------------------------------
# Same call, the apoptosis gene set. Adds a column "ApopScore1".
# Hint: AddModuleScore(); args features = list(apoptosis_genes), name = "ApopScore"
# YOUR CODE HERE:
combined <-

# Compare the suspect sample to everyone else. (written for you)
p_stress <- VlnPlot(combined, features = c("StressScore1", "ApopScore1", "percent.mt"),
                    group.by = "sample", pt.size = 0) +
  patchwork::plot_annotation(title = "Stress / apoptosis / mito by sample")
ggsave(file.path(OUT_DIR, "gata1_module_scores_by_sample.png"),
       plot = p_stress, width = 12, height = 5, dpi = 150)

# ---- 7. Save the filtered, annotated object for script 02 -------------------
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_qc.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_qc.rds"))
