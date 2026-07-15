# =============================================================================
# 03_cell_type_annotation.R
# -----------------------------------------------------------------------------
# Step 3: automatically annotate cell types. Three complementary approaches:
#   A. Reference mapping with Seurat (transfer labels from a labeled reference)
#   B. Reference-based scoring with SingleR + celldex (built-in atlases)
#   C. Marker-gene scoring with SC-Type (no labeled reference needed)
#
# DATASET: the pre-labeled T21 object (Trisomy 21 fetal bone marrow). Because it
# already has TRUE labels (cell.labels / broad_extfig7A_cell.labels) we can split
# it in half and measure how well each method recovers the known answer.
#
# Companion lesson: lessons/03_cell_type_annotation.md
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(ggplot2)
library(dplyr)
library(HGNChelper)   # used by SC-Type to fix gene symbols
library(openxlsx)     # used by SC-Type to read its marker database
library(SingleR)      # reference-based annotation
library(celldex)      # ready-made reference atlases for SingleR

# -----------------------------------------------------------------------------
# 0. Load the labeled object
# -----------------------------------------------------------------------------
t21 <- load_one(T21_LABELED)   # robustly grabs whatever object the file holds
ls()

# The reduction holding the published UMAP may be called "Xumap_" or "umap";
# detect it so the plots below work either way.
umap_name <- if ("Xumap_" %in% Reductions(t21)) "Xumap_" else "umap"
# The fine-grained label column is "cell.labels"; the broad one is below.
broad_col <- if ("broad_extfig7A_cell.labels" %in% colnames(t21@meta.data))
  "broad_extfig7A_cell.labels" else "cell.labels"

# Split the labeled data in half: a "reference" we pretend is annotated and a
# "query" we pretend is unknown. Seed keeps everyone's split identical.
set.seed(42)
query_cells <- c()
for (cell_state in unique(t21@meta.data[[broad_col]])) {
  filt <- t21@meta.data[t21@meta.data[[broad_col]] == cell_state, ]
  query_cells <- c(query_cells, sample(rownames(filt), nrow(filt) / 2))
}
ref_cells <- setdiff(rownames(t21@meta.data), query_cells)

t21_query <- subset(t21, cells = query_cells)
t21_ref   <- subset(t21, cells = ref_cells)
rm(t21)

# =============================================================================
# A. SEURAT REFERENCE MAPPING
# =============================================================================
# 1. FindTransferAnchors: pairs of mutual-nearest-neighbor cells between the
#    reference and query that anchor the label transfer.
t21_query.anchors <- FindTransferAnchors(reference = t21_ref, query = t21_query)

# 2. TransferData: use the anchors to predict a label for every query cell.
ls_predictions <- TransferData(anchorset = t21_query.anchors,
                               refdata = t21_ref@meta.data[[broad_col]])

# 3. AddMetaData: attach the predictions to the query object.
identical(rownames(t21_query@meta.data), rownames(ls_predictions))  # sanity check
t21_query@meta.data$ls_predicted_id <- ls_predictions$predicted.id

# 4. Score it: since we know the truth, count how often the prediction matches.
t21_query$ls_prediction.match <- t21_query$ls_predicted_id == t21_query@meta.data[[broad_col]]
table(t21_query$ls_prediction.match)

# =============================================================================
# B. SingleR + celldex (reference-based, using a built-in atlas)
# =============================================================================
# celldex provides curated reference datasets. For bone marrow / blood, the
# Human Primary Cell Atlas or Monaco immune reference are good choices.
hpca <- celldex::HumanPrimaryCellAtlasData()

# SingleR correlates each query cell's expression against the reference's per-
# label profiles and assigns the best match. Feed it the log-normalized data.
singleR_pred <- SingleR(
  test   = GetAssayData(t21_query, assay = "RNA", layer = "data"),
  ref    = hpca,
  labels = hpca$label.main
)
t21_query$singleR_id <- singleR_pred$labels
table(t21_query$singleR_id)

# =============================================================================
# C. SC-TYPE (marker-gene scoring, no labeled reference required)
# =============================================================================
# 1. Load SC-Type's functions and its marker database directly from GitHub.
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_    <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Immune system"                 # pick the tissue closest to your data
gs_list <- gene_sets_prepare(db_, tissue) # positive + negative marker sets
names(gs_list$gs_positive)                # the cell types this DB can call

# 2. SC-Type scores cells on SCALED data. Extract the scaled matrix (v4/v5 safe).
seurat_v5 <- isFALSE("counts" %in% names(attributes(t21_query[["RNA"]])))
scaled <- if (seurat_v5) as.matrix(t21_query[["RNA"]]$scale.data)
          else           as.matrix(t21_query[["RNA"]]@scale.data)

# 3. Score every cell against every candidate cell type.
es.max <- sctype_score(scRNAseqData = scaled, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# 4. Aggregate scores per CLUSTER and keep the top-scoring type for each.
#    Here we use the known broad labels as "clusters"; in a blind run you would
#    use seurat_clusters instead.
cL_results <- do.call("rbind", lapply(unique(t21_query@meta.data[[broad_col]]), function(cl) {
  es.max.cl <- sort(rowSums(es.max[, rownames(t21_query@meta.data[t21_query@meta.data[[broad_col]] == cl, ])]),
                    decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                  ncells = sum(t21_query@meta.data[[broad_col]] == cl)), 10)
}))
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# Flag low-confidence calls (score small relative to cell count).
sctype_scores$confidence <- "High"
sctype_scores$confidence[as.numeric(as.character(sctype_scores$scores)) <
                           sctype_scores$ncells / 4] <- "Low"
sctype_scores

# 5. Write the SC-Type label back onto each cell and plot.
t21_query@meta.data$sctype_classification <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  t21_query@meta.data$sctype_classification[t21_query@meta.data[[broad_col]] == j] <-
    as.character(cl_type$type[1])
}

DimPlot(t21_query, reduction = umap_name, group.by = "sctype_classification",
        label = TRUE, repel = TRUE)
DimPlot(t21_query, reduction = umap_name, group.by = broad_col,
        label = TRUE, repel = TRUE)
