source("~/srworkshop/projectA/00_paths_and_setup.R")

library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

cellChat <- readRDS(file.path(OUT_DIR, "gata1_cellchat.rds"))
#If you didn't make this yet, use mine!
#cellChat <- readRDS(file.path(COOKING, "gata1_cellchat.rds"))



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
# mergeCellChat() + comparison plots. This is a optional exersize. 
