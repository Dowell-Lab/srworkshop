# =============================================================================
# GATA1 track - Script 03: Multi-reference SingleR annotation + marker genes
# -----------------------------------------------------------------------------
# Assigns a cell-type label to every cell automatically with SingleR, and
# shows WHY you might consult more than one reference atlas for hematopoiesis.
# Then sanity-checks the labels against canonical marker genes.
#
# Inputs : OUT_DIR/gata1_combined_clustered.rds  (from script 02)
# Outputs: an annotated object + UMAPs / FeaturePlots in OUT_DIR.
# -----------------------------------------------------------------------------
# HOW TO USE THIS TEMPLATE
#   Work through it top to bottom alongside the worksheet
#   (03_multiref_singler_annotation.md). The boilerplate is already written for
#   you. Wherever you see a block like:
#
#       # ---- Step 1a: one matrix for SingleR ... ----
#       # Hint: JoinLayers(); args ...
#       # YOUR CODE HERE:
#       combined_joined <-
#
#   finish the line yourself. If you get stuck, the completed answer key is in
#   scripts_finished/03_multiref_singler_annotation.R -- try it on your own first.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(patchwork)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_clustered.rds"))
#If you didn't make this yet, use mine!
#combined <- readRDS(file.path(COOKING, "gata1_combined_clustered.rds"))


# ---- 1. Prepare a single expression matrix for SingleR ----------------------
# After merge(), the RNA assay holds one "data" layer PER sample
# (data.EuploidGATA1sD7, data.T21wtGATA1D9, ...). SingleR wants ONE matrix,
# so JoinLayers stitches them into a unified "data" layer first.
Layers(combined[["RNA"]])           # inspect the per-sample layers

# ---- Step 1a: join the per-sample layers ------------------------------------
# Stitch the per-sample "data" layers into one unified layer so SingleR sees
# every cell in a single matrix.
# Hint: JoinLayers(); args combined
# YOUR CODE HERE:
combined_joined <-

# ---- Step 1b: pull out the log-normalized matrix ----------------------------
# Grab the joined "data" layer (log-normalized values) -- that is exactly what
# SingleR's correlation scoring expects.
# Hint: LayerData(); args combined_joined, assay = "RNA", layer = "data"
# YOUR CODE HERE:
sce_counts <-

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




# This tiny wrapper is written for you: it correlates every cell against the
# reference's log-counts and returns the best-matching label per cell.
run_singler <- function(ref, test_mat) {
  SingleR(
    test   = test_mat,
    ref    = assay(ref, "logcounts"),
    labels = ref$label.main
  )
}

# ---- Step 2a: run SingleR against the broad reference -----------------------
# Score every cell in your matrix against ref_broad using the helper above.
# Hint: run_singler(); args ref_broad, sce_counts
# YOUR CODE HERE:
pred_broad <-
table(pred_broad$labels)

# ---- 3. Collapse rare labels into "other" -----------------------------------
# With many references you get a long tail of labels assigned to only a handful
# of cells. Those clutter the UMAP. We keep labels with > 150 cells and bucket
# the rest as "other" for a readable plot, while keeping the full label too.

# ---- Step 3a: bucket rare labels --------------------------------------------
# add_count() counts how many cells share each label; then case_when keeps the
# label when it is common (> 150 cells) and relabels the rest as "other".
# Hint: mutate(); labels_other = case_when(label_n > 150 ~ labels, TRUE ~ "other")
# YOUR CODE HERE:
pred_df <- as.data.frame(pred_broad) %>%
  add_count(labels, name = "label_n") %>%
  mutate(labels_other = )

# Write the label versions back into the Seurat metadata. Three of these are
# done for you; you add the "_other" column that CellChat (script 04) will use.
combined$SingleR_label        <- pred_df$labels
# ---- Step 3b: the readable-label column (script 04 depends on this) ---------
# Store the bucketed labels so the next lesson can group cells by cell type.
# Hint: assign pred_df$labels_other into combined$SingleR_labels_other
# YOUR CODE HERE:
combined$SingleR_labels_other <-
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

# ---- 5. Visualize the automated annotation ----------------------------------
save_dim <- function(p, file, w = 8, h = 6) {
  ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
}

# ---- Step 5a: UMAP colored by the full SingleR label (WORKED EXAMPLE) -------
# This first DimPlot is done for you. Notice the pattern:
#   DimPlot(object, reduction = "umap", group.by = <a metadata column>)
# then print and save it. Copy this pattern for the two below, changing only
# group.by.
p<- DimPlot(combined, reduction = "umap", group.by = "SingleR_label", label = TRUE)
p
fn = paste0(whichlabelref, "gata1_umap_singler_label.png")
save_dim(p,
         fn)

# ---- Step 5b: UMAP colored by the pruned (high-confidence) label ------------
# Hint: same DimPlot pattern, group.by = "SingleR_pruned"
# YOUR CODE HERE:
p<-
p
fn <- paste0(whichlabelref, "gata1_umap_singler_pruned.png")

save_dim(p,
         fn)

# ---- Step 5c: UMAP colored by the readable "_other" label -------------------
# Hint: group.by = "SingleR_labels_other"
# YOUR CODE HERE:
p<-
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

# One marker at a time is done for you as the worked example: FeaturePlot paints
# a gene's expression straight onto the UMAP.
agene = "HBG2"

p_markers <- FeaturePlot(combined, features = agene, reduction = "umap",
                         pt.size = 0.3, ncol = 3)
p_markers
fn = paste0("gata1dataset_", agene,".png")
save_dim(p_markers, fn, w = 12, h = 14)


markers <- c("HBG2", "HBG1", "GATA1", "GFI1B", "TFRC",
             "CD34", "PTPRC", "ITGA4", "VAMP8", "DIAPH3")

# ---- Step 6a: FeaturePlot the whole marker panel ----------------------------
# Same FeaturePlot call as above, but hand it the full `markers` vector so every
# gene gets its own UMAP panel.
# Hint: FeaturePlot(); args combined, features = markers, reduction = "umap",
#       pt.size = 0.3, ncol = 3
# YOUR CODE HERE:
p_markers <-
p_markers
save_dim(p_markers, "gata1_feature_markers.png", w = 12, h = 14)

# Side-by-side: labels next to a key marker. (written for you)
p_combo <- DimPlot(combined, reduction = "umap", group.by = "SingleR_label", label = TRUE) +
  FeaturePlot(combined, features = "HBG1", reduction = "umap")
p_combo
save_dim(p_combo, "gata1_labels_vs_HBG1.png", w = 13, h = 6)

# ---- 7. Save the fully annotated object -------------------------------------
# These saves are commented out. Uncomment them to write the annotated objects.
# script 04 (CellChat) reads gata1_combined_annotated_joined.rds, so uncomment
# at least the second pair before moving on.
#saveRDS(combined, file.path(OUT_DIR, "gata1_combined_annotated.rds"))
#message("Saved: ", file.path(OUT_DIR, "gata1_combined_annotated.rds"))

#saveRDS(combined_joined, file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
#message("Saved: ", file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))

rm(combined_joined)
