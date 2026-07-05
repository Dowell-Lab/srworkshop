# =============================================================================
# GATA1 track - Script 04: Cell-cell communication with CellChat
# -----------------------------------------------------------------------------
# Infers signaling between the cell types you annotated in Lesson 03, using the
# GATA1 (GSE271399) object you built yourself. Unlike the fetal-bone-marrow
# track (which loads a PRE-MADE CellChat .rds), here you build the CellChat
# object end to end from your own annotated Seurat object.
#
# Input : OUT_DIR/gata1_combined_annotated.rds  (from script 03)
# Output: a CellChat object + signaling plots in OUT_DIR.
#
# Based on the CellChat vignette:
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))

combined #look at how many cells you have now

gc() # garbage collector, cleans up how much memory the script is using to a minimum 

# ---- 1. Choose the grouping label -------------------------------------------
# CellChat groups cells by a cell-type label. We use the SingleR labels from
# script 03. Use the "_other" version so rare labels don't create tiny groups
# that CellChat can't model. Drop cells with no usable label.

combined$cc_group <- combined$SingleR_labels_other
combined <- subset(combined, subset = !is.na(cc_group) & cc_group != "other")

combined #look at how many cells you have now

data.input <- LayerData(
  combined,
  assay = "RNA",
  layer = "data"   # or "scale.data"/"counts" depending on your workflow
)  # returns a matrix‑like object[][]

# Build meta with your cc_group labels
meta <- data.frame(
  cc_group = combined$cc_group,
  row.names = colnames(combined)
)

# ---- 2. Build the CellChat object -------------------------------------------

cellChat <- createCellChat(
  object   = data.input,
  meta     = meta,
  group.by = "cc_group"
)

# Attach the human ligand-receptor database.
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellChat@DB <- CellChatDB

# ---- 3. Pre-processing: over-expressed genes & interactions -----------------
cellChat <- subsetData(cellChat)            # restrict to genes in the DB
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

# ---- 4. Infer the communication network (the slow step) ---------------------
# computeCommunProb estimates the probability of signaling between each pair of
# cell groups for every ligand-receptor pair. On a small server this can take a
# while; run once, then save the object so you never recompute.
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)  # drop tiny groups
cellChat <- computeCommunProbPathway(cellChat)             # roll up to pathways
cellChat <- aggregateNet(cellChat)                         # count edges

saveRDS(cellChat, file.path(OUT_DIR, "gata1_cellchat.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_cellchat.rds"))

# ---- 5. Whole-dataset overview ----------------------------------------------
# Circle plot: nodes = cell types, edge thickness = communication strength.
groupSize <- as.numeric(table(cellChat@idents))

png(file.path(OUT_DIR, "gata1_cellchat_overview_circle.png"),
    width = 1400, height = 700, res = 150)
par(mfrow = c(1, 2))
netVisual_circle(cellChat@net$count,  vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction strength")
dev.off()

# ---- 6. Which pathways drive signaling? -------------------------------------
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming")

png(file.path(OUT_DIR, "gata1_cellchat_signaling_roles.png"),
    width = 1400, height = 800, res = 150)
ComplexHeatmap::draw(ht1 + ht2)
dev.off()

# ---- 7. Focus on one pathway ------------------------------------------------
# Pick a pathway that is actually present in YOUR result (the set depends on the
# cell types recovered). List what's available, then visualize one.
available_pathways <- cellChat@netP$pathways
print(available_pathways)

if (length(available_pathways) > 0) {
  pathways.show <- available_pathways[1]   # swap for a pathway of interest

  png(file.path(OUT_DIR, "gata1_cellchat_pathway_circle.png"),
      width = 800, height = 800, res = 150)
  netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")
  dev.off()

  png(file.path(OUT_DIR, "gata1_cellchat_pathway_heatmap.png"),
      width = 800, height = 700, res = 150)
  netVisual_heatmap(cellChat, signaling = pathways.show, color.heatmap = "Reds")
  dev.off()

  # Contribution of each ligand-receptor pair within the pathway.
  p_contrib <- netAnalysis_contribution(cellChat, signaling = pathways.show)
  ggplot2::ggsave(file.path(OUT_DIR, "gata1_cellchat_pathway_contribution.png"),
                  plot = p_contrib, width = 7, height = 4, dpi = 150)
}

# ---- 8. (Optional) compare conditions ---------------------------------------
# To compare, say, GATA1s vs wtGATA1: subset the annotated object by construct,
# build a CellChat object for each (steps 2-4), then use CellChat's
# mergeCellChat() + comparison plots. Left as an exercise.
