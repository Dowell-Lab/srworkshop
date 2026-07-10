# =============================================================================
# 00_paths_and_setup.R
# -----------------------------------------------------------------------------
# Central place for the file paths on the AWS / Fiji cluster.
# The data files are READ-ONLY and shared, so nobody downloads anything.
#
# `source("00_paths_and_setup.R")` at the top of any lesson script to get
# these paths. If the workshop directory ever moves, edit ONLY this file.
#
# This course has TWO tracks:
#   - GATA1 (PRIMARY) : GSE271399 iPSC -> hematopoiesis time course, run end to
#                       end by students (scripts/gata1/).
#   - Fetal bone marrow (REFERENCE) : T21 (Trisomy 21) and D21 (disomic control)
#                       samples for a few extra techniques (scripts/fbm/).
# =============================================================================

# ---- Root of the shared workshop data ----
SREAD_ROOT <- "/scratch/Shares/public/sread/cookingShow"

# =============================================================================
# PRIMARY DATASET: GATA1 / GSE271399 (iPSC -> hematopoiesis time course)
# -----------------------------------------------------------------------------
# Takasaki et al., Stem Cell Reports 2025: "Single-cell transcriptomics reveal
# individual and cooperative effects of trisomy 21 and GATA1s on hematopoiesis."
# Isogenic euploid vs T21 iPSC lines, each carrying wild-type GATA1 or the
# truncated GATA1s isoform, sampled at differentiation days D7/D9/D11.
# This is the dataset students run END TO END.
# =============================================================================
GATA1_DIR <- file.path(SREAD_ROOT, "day7a", "iPCStoblood")

# The 12 condition names. Files follow GSE271399_<sample>_{matrix.mtx,barcodes.tsv,features.tsv}.gz
GATA1_SAMPLES <- c(
  "EuploidGATA1sD7",  "EuploidGATA1sD9",  "EuploidGATA1sD11",
  "EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
  "T21GATA1sD7",      "T21GATA1sD9",      "T21GATA1sD11",
  "T21wtGATA1D7",     "T21wtGATA1D9",     "T21wtGATA1D11"
)


#GATA1_SAMPLES <- c(
#  "EuploidGATA1sD7",  "EuploidGATA1sD11",
#  "EuploidwtGATA1D7", "EuploidwtGATA1D11",
#  "T21GATA1sD7",      "T21GATA1sD11",
#  "T21wtGATA1D7",     "T21wtGATA1D11"
#)

# =============================================================================
# REFERENCE DATASET: T21 / D21 fetal bone marrow (the original demo track)
# =============================================================================

# ---- 1. Raw 10x CellRanger matrices (QC / integration chapter) ----
# Each sample folder contains outs/filtered_feature_bc_matrix/
RAW10X_DIR   <- file.path(SREAD_ROOT, "day7a", "total_data")
T21_SAMPLE   <- "T21BM_male04"   # the Trisomy 21 sample we standardize on
D21_SAMPLE   <- "D21_male35"     # the Disomic control sample we standardize on

# ---- 2. Pre-labeled, normalized Seurat objects ----
# Used for annotation, pseudotime, and CellChat Part 1.
# These carry the published cell-type labels and a UMAP embedding.
LABELED_DIR  <- file.path(SREAD_ROOT, "day9", "cellchat-prelabeled-seurat-objs")
T21_LABELED  <- file.path(LABELED_DIR, "t21_norm_seurat_obj.RData")
D21_LABELED  <- file.path(LABELED_DIR, "d21_norm_seurat_obj.RData")

# ---- 3. Pre-made CellChat objects (CellChat Part 2) ----
# These are .rds files, so load them with readRDS().
CELLCHAT_DIR    <- file.path(SREAD_ROOT, "day9", "cellchatobjs")
CELLCHAT_T21_RDS <- file.path(CELLCHAT_DIR, "t21-cellchat-obj.rds")
CELLCHAT_D21_RDS <- file.path(CELLCHAT_DIR, "d21-cellchat-obj.rds")

# ---- Where YOU can write intermediate results (NOT the read-only share) ----
# Defaults to a folder in your home directory; change as you like.
username = Sys.info()[["user"]] 



OUT_DIR <- path.expand(paste0("/localscratch/Users/",username,"/scRNAseq_course_output"))
COOKING <- path.expand(paste0("/scratch/Shares/public/sread/cookingShow/singlecellgeneral/"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Small helper: load an .Rdata file and return the object it contains ----
# The pre-labeled objects may load under a name like "t21" or "seurat_obj";
# this grabs whatever single object the file holds so you don't have to guess.
load_one <- function(path) {
  e <- new.env()
  nm <- load(path, envir = e)        # nm = names of objects that were loaded
  if (length(nm) == 1) return(get(nm, envir = e))
  # If more than one object, print the names so you can pick.
  message("Multiple objects in ", basename(path), ": ", paste(nm, collapse = ", "))
  mget(nm, envir = e)
}

