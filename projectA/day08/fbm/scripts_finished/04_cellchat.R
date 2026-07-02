# =============================================================================
# 04_cellchat.R
# -----------------------------------------------------------------------------
# Step 4: infer cell-cell communication with CellChat.
#   PART 1 : build a CellChat object from a Seurat object (slow steps shown).
#   PART 2 : analyze a PRE-MADE CellChat object (fast; what we run in class).
#
# DATASET: the T21 (Trisomy 21) fetal bone marrow object. A matching D21 control
# object exists at the same paths if you want to compare conditions afterward.
#
# Based on the CellChat vignette:
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
#
# Companion lesson: lessons/04_cellchat.md
# =============================================================================

source("../00_paths_and_setup.R")   # gives T21_LABELED, CELLCHAT_T21_RDS, load_one()

library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

# =============================================================================
# PART 1 — Building a CellChat object (slow; usually done ahead of time)
# =============================================================================
# Run this section once on your own time. In class, skip to PART 2 and load the
# pre-made object instead.

t21 <- load_one(T21_LABELED)

## 2. Create a CellChat object, grouping cells by their cell-type labels.
cellChat <- createCellChat(object = t21, group.by = "cell.labels")

## 3. Attach the human ligand-receptor database.
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)          # see what's in the DB
dplyr::glimpse(CellChatDB$interaction)    # peek at its structure

## 4. (Optional) Restrict to "Secreted Signaling" to save memory.
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellChat@DB <- CellChatDB        # using the full DB here; swap in .use to subset
cellChat <- subsetData(cellChat)

## 5. Identify over-expressed genes and ligand-receptor pairs (a few minutes).
future::plan("multisession", workers = 4)   # parallelize to ease memory
cellChat <- identifyOverExpressedGenes(cellChat, do.fast = FALSE)
cellChat <- identifyOverExpressedInteractions(cellChat)

## 6. Infer the communication network (VERY slow on a small server).
# Left commented because we use the pre-made object below in class.
# cellChat <- computeCommunProb(cellChat, type = "triMean")

# =============================================================================
# PART 2 — Analyzing a PRE-MADE CellChat object (fast; run this in class)
# =============================================================================
# Start fresh here: load the pre-computed object (.rds -> readRDS).
cellChat <- readRDS(CELLCHAT_T21_RDS)

## 2. Drop communications supported by too few cells (min 10 per group).
cellChat <- filterCommunication(cellChat, min.cells = 10)

## 3. Aggregate single L-R communication up to the pathway level.
cellChat <- computeCommunProbPathway(cellChat)

## 4. Count edges (interactions) between every pair of cell types.
cellChat <- aggregateNet(cellChat)

## 5. Whole-dataset overview: a circle plot of all interactions.
# Nodes = cell types; edge thickness = communication strength.
groupSize <- as.numeric(table(cellChat@idents))
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction weights/strength")

## 6. Which pathways matter? Network-centrality heatmaps.
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming")
ht1 + ht2   # darker = stronger signaling for that pathway x cell type

## 7. Focus on one pathway. In Trisomy 21, interferon-related signaling is of
##    special interest (4 of the ~225 chr-21 genes are IFN-related). We look at
##    CCL as an example here.
pathways.show <- c("CCL")

# a. Circle plot for just this pathway.
par(mfrow = c(1, 1))
netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")

# b. Heatmap for this pathway.
par(mfrow = c(1, 1))
netVisual_heatmap(cellChat, signaling = pathways.show, color.heatmap = "Reds")

# c. Contribution of each ligand-receptor pair within the pathway.
netAnalysis_contribution(cellChat, signaling = pathways.show)

# d. Bubble plot of significant L-R pairs from chosen senders to receivers.
levels(cellChat@idents)   # see the index for each cell type
netVisual_bubble(cellChat, sources.use = 4, targets.use = c(5:7), remove.isolate = FALSE)

# Limit the bubble plot to a few pathways of interest.
netVisual_bubble(cellChat, sources.use = 4, targets.use = c(5:7),
                 signaling = c("CCL", "CXCL", "TNF"), remove.isolate = FALSE)

# To COMPARE conditions later, repeat PART 2 with CELLCHAT_D21_RDS and use
# CellChat's mergeCellChat() + comparison plots.
