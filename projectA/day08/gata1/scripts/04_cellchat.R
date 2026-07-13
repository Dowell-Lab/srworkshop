# =============================================================================
# GATA1 track - Script 04: Cell-cell communication with CellChat
# -----------------------------------------------------------------------------
# Infers signaling between the cell types you annotated in Lesson 03, using the
# GATA1 (GSE271399) object you built yourself. Unlike the fetal-bone-marrow
# track (which loads a PRE-MADE CellChat .rds), here you build the CellChat
# object end to end from your own annotated Seurat object.
#
# Input : OUT_DIR/gata1_combined_annotated_joined.rds  (from script 03)
# Output: a CellChat object + signaling plots in OUT_DIR.
#
# Based on the CellChat vignette:
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
# -----------------------------------------------------------------------------
# HOW TO USE THIS TEMPLATE
#   Work through it top to bottom alongside the worksheet (04_cellchat.md). The
#   boilerplate is already written for you. Wherever you see a block like:
#
#       # ---- Step 2a: build the CellChat object ... ----
#       # Hint: createCellChat(); args ...
#       # YOUR CODE HERE:
#       cellChat <-
#
#   finish the line yourself. If you get stuck, the completed answer key is in
#   scripts_finished/04_cellchat.R -- try it on your own first.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
#If you didn't make this yet, use mine!
#combined <- readRDS(file.path(COOKING, "gata1_combined_annotated_joined.rds"))


combined #look at how many cells you have now

gc() # garbage collector, cleans up how much memory the script is using to a minimum

# ---- 1. Choose the grouping label -------------------------------------------
# CellChat groups cells by a cell-type label. We use the SingleR labels from
# script 03. Use the "_other" version so rare labels don't create tiny groups
# that CellChat can't model. Drop cells with no usable label.

# ---- Step 1a: pick the grouping label ---------------------------------------
# Copy the readable SingleR labels into a new "cc_group" column CellChat will
# group by.
# Hint: assign combined$SingleR_labels_other into combined$cc_group
# YOUR CODE HERE:
combined$cc_group <-

# ---- Step 1b: drop the unusable cells ---------------------------------------
# Keep only cells that got a real label -- no NAs, and not the "other" bucket --
# so CellChat isn't asked to model groups too small to be meaningful.
# Hint: subset(); args combined, subset = !is.na(cc_group) & cc_group != "other"
# YOUR CODE HERE:
combined <-

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

# ---- Step 2a: create the CellChat object ------------------------------------
# Hand CellChat the log-normalized matrix, the meta table, and the column to
# group cells by.
# Hint: createCellChat(); args object = data.input, meta = meta,
#       group.by = "cc_group"
# YOUR CODE HERE:
cellChat <-

# Attach the human ligand-receptor database.
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# ---- Step 2b: attach the ligand-receptor database ---------------------------
# Point the CellChat object at the human L-R database so it knows which
# ligand-receptor pairs to score.
# Hint: assign CellChatDB into the object's @DB slot: cellChat@DB
# YOUR CODE HERE:
cellChat@DB <-

# ---- 3. Pre-processing: over-expressed genes & interactions -----------------
# ---- Step 3a: restrict to genes in the database -----------------------------
# Drop every gene that isn't part of a signaling interaction.
# Hint: subsetData(); args cellChat
# YOUR CODE HERE:
cellChat <-

# ---- Step 3b: find genes over-expressed in each group -----------------------
# Hint: identifyOverExpressedGenes(); args cellChat
# YOUR CODE HERE:
cellChat <-

# ---- Step 3c: keep L-R pairs where both partners are over-expressed ---------
# Hint: identifyOverExpressedInteractions(); args cellChat
# YOUR CODE HERE:
cellChat <-

# ---- 4. Infer the communication network (the slow step) ---------------------
# computeCommunProb estimates the probability of signaling between each pair of
# cell groups for every ligand-receptor pair. On a small server this can take a
# while; run once, then save the object so you never recompute.

# ---- Step 4a: compute communication probabilities (SLOW) --------------------
# HEADS UP: this is the expensive call -- it can take a while. Run it once.
# Hint: computeCommunProb(); args cellChat, type = "triMean"
# YOUR CODE HERE:
cellChat <-

# ---- Step 4b: drop edges from tiny groups -----------------------------------
# Hint: filterCommunication(); args cellChat, min.cells = 10
# YOUR CODE HERE:
cellChat <-

# ---- Step 4c: roll individual L-R pairs up to pathway level -----------------
# Hint: computeCommunProbPathway(); args cellChat
# YOUR CODE HERE:
cellChat <-

# ---- Step 4d: count edges to build the summary network ----------------------
# Hint: aggregateNet(); args cellChat
# YOUR CODE HERE:
cellChat <-

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
# ---- Step 7a: list the pathways your run recovered --------------------------
# Pull the available pathways out of the object so we never hard-code one that
# might be absent.
# Hint: read the netP slot: cellChat@netP$pathways
# YOUR CODE HERE:
available_pathways <-
print(available_pathways)

# The per-pathway loop is written for you: it walks every recovered pathway and
# saves a circle plot + an L-R contribution bar chart for each.
if (length(available_pathways) > 0) {for (i in seq(length(available_pathways))){
  onepathway = available_pathways[i]
  pathways.show <- available_pathways[i]   # swap for a pathway of interest
  fn = paste0(onepathway, "_gata1_cellchat_pathway_circle.png")
  png(file.path(OUT_DIR, fn),
      width = 800, height = 800, res = 150)
  netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")
  dev.off()

  # Contribution of each ligand-receptor pair within the pathway.
  p_contrib <- netAnalysis_contribution(cellChat, signaling = pathways.show)
  fn = paste0(onepathway, "_gata1_cellchat_pathway_contribution.png")
  ggplot2::ggsave(file.path(OUT_DIR, fn),
                  plot = p_contrib, width = 7, height = 4, dpi = 150)
}}

# ---- 8. (Optional) compare conditions ---------------------------------------
# To compare, say, GATA1s vs wtGATA1: subset the annotated object by construct,
# build a CellChat object for each (steps 2-4), then use CellChat's
# mergeCellChat() + comparison plots. This is a optional exersize.
