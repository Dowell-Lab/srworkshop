# Author: Georgia Barone, 2024
# Day 08, CellChat Part 2 Analysis with scRNA-seq data
# Note: This worksheet is based off of the CellChat vignette: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

# Load in required libraries 
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

## 1. Load in pre-made CellChat object. 
# Note** When performing your own analysis you will need to make your own CellChat object (directions on how to do so are in Cell-Chat-Part-1.R)
# For today you can use the object provided below

# Load in the pre-made trisomic CellChat object
load("/scratch/Shares/public/sread2024/cookingShow/day9/cellchat/cellchat_object-t21.RData")

# If you have completed the worksheet and are running the disomic object pull in this file instead: /scratch/Shares/public/sread2024/cookingShow/day9/cellchat/cellchat_object-d21.RData

## 2. Set the minimum number of cells required to assess cell-cell communication on a per-network basis
# We will set the min number of cells required in each cell group to be 10
cellChat <- filterCommunication(cellChat, min.cells = 10)

## 3. Recalculate cell-cell communication (after filtering out networks with less than 10 interacting cell types)
cellChat <- computeCommunProbPathway(cellChat)

## 4. Count the number of edges (cell-cell communications) on a per-network basis
cellChat <- aggregateNet(cellChat)

## 5. Visualize ALL cell-cell communication within a dataset
# While CellChat is best used for analysis of cell-cell communication within a specific pathway, you are also able to visualize cell-cell communication across an entire dataset
# The code below will output a circle plot, where the nodes are cell types and edges represent the communications between cell types
# The thicker an edge is, the more likely two cell types are to be interacting with each other (higher communication probability between two cell types)

# Set the group size (below we have set the group to represent all cell identities)
groupSize <- as.numeric(table(cellChat@idents))

# Plot the circle plot
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# As you can see, there is not a ton of information we can pull from this graph, because all interactions are being plotted on top of each other
# But how do we know which pathway to look at?

# Answer: 
## 6. Let's make a heatmap of signals contributing the most to outgoing or incoming signaling of certain cell groups
# First let's compute network centrality scores (see slides)
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP") # 'netP' means the inferred intercellular communication network of signaling pathways

# create outgoing heat map
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing")

# create outgoing heat map
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming")

# show both heatmaps next to each other
ht1 + ht2

# This plot helps to identify which pathways you should analyze. 
# How dark the green in the heat map indicates how strong the outgoing or incoming signal is for a certain pathway (y axis) and certain cell type (x axis)
# This is a good tool to utilize for future exploratory analysis! If you are unsure which pathways to look at in real life - start with these heat maps

# 7. Plot on a per-pathway basis
# 4 of the 225 genes on chromosome 21 are interferon (IFN)-related genes [1]. 
# Individuals with Trisomy 21 have an extra 21st chromosome, making IFN-related pathways a point of interest.

# Below are some IFN-related pathways.
# ifn_paths: "IFN-II", "MHC-I", "MHC-II", "CXCL", "CCL", "CD40", "TNF"

# Let's examine the CCL pathway in a variety of different ways
# Setting the pathway to 'CCL'
pathways.show <- c("CCL")

# a. Let's start with a circle plot, but this time specific to the CCL pathway
par(mfrow=c(1,1)) # set the plotting layout to a single plot
netVisual_aggregate(cellChat, signaling = pathways.show, layout = "circle")

# b. Let's plot a heatmap now for the CCL pathway
par(mfrow=c(1,1)) # set the plotting layout to a single plot
netVisual_heatmap(cellChat, signaling = pathways.show, color.heatmap = "Reds")

# c. Let's compute the contribution of each ligand-receptor pair to the overall signaling pathway 
# We can visualize the contributions of each L-R within a certain signaling pathway
netAnalysis_contribution(cellChat, signaling = pathways.show)

# d. Plot all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# To know what cell types you have selected, run:
levels(cellChat@idents)

# Here, I am setting "CD56 bright NK" [4] as the senders; and "CMP [5] "DC1" [6] "DC2"[7] as the receivers (see x-axis)
netVisual_bubble(cellChat, sources.use = 4, targets.use = c(5:7), remove.isolate = FALSE)

# We can also limit which pathways are shown in the bubble plot
# Here we will create the same plot as above (showing significant L-R pairs), but limit the number of pathways shown ("CCL","CXCL","TNF")
netVisual_bubble(cellChat, sources.use = 4, targets.use = c(5:7), signaling = c("CCL","CXCL","TNF"), remove.isolate = FALSE)

# Nice work! Feel free to make more plots by following the CellChat vignette here: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

# Some suggestions for plots (there are tutorials for each of these in the CellChat vignette):
# 1. Plot the signaling gene expression distribution using violin/dot plot 
# 2. Visualize the computed centrality scores using a heat map (similar to step 12, but for a specific pathway)
# 3. Visualize dominant senders (sources) and receivers (targets) in a 2D space
# 4. Identify and visualize outgoing communication pattern of secreting cells

# Refs:
# 1. Chung, H., Green, P. H., Wang, T. C., & Kong, X. F. (2021). Interferon-driven immune dysregulation in Down syndrome: a review of the evidence. Journal of Inflammation Research, 5187-5200.
# 2. https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html



