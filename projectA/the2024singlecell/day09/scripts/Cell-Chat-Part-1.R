# Author: Georgia Barone, 2024
# Day 08, CellChat Part 1 Analysis with scRNA-seq data
# Note: This worksheet is based off of the CellChat vignette: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

# Load in required libraries 
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

## 1. Load in pre-made Seurat object. 
# Note** When performing your own analysis you will need to make your own Seurat object, but for today you can use the object provided below

# Load in the Seurat object - let's start with the Trisomic sample (t21_official_umap_clust_06.25.24.Rdata)
# Note** There is also a Disomic sample in this directory you are welcome to run CellChat on for comparison after completing both worksheets (/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/d21_subset_official_umap_clust_06.25.24.Rdata))

load("/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/t21_official_umap_clust_06.25.24.Rdata")

## 2. Create CellChat object (from Seurat object)
# Note: a warning here about the 'samples' column is okay, proceed
cellChat <- createCellChat(object = t21, group.by = "cell.labels")

## 3. Initialize CellChat database
CellChatDB <- CellChatDB.human

# This line shows you what makes up the CellChat database
showDatabaseCategory(CellChatDB)

# Take a look at the CellChat database structure
dplyr::glimpse(CellChatDB$interaction)

## 4. Add CellChatDB to the CellChat object
# We will be using the Secreted Signaling subset of the CellChat DB to save memory on the AWS today
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellChat@DB <- CellChatDB
cellChat <- subsetData(cellChat)

## 5. Identify highly variable ligand/receptor pairs
# In this step, CellChat identifies highly variable ligand/receptor pairs (differentially expressed across prelabeled celltypes)
future::plan("multisession", workers = 4) # run this in parallel, helps with memory allocation 

# Note ** The following line will take 3 - 5 mins
cellChat <- identifyOverExpressedGenes(cellChat, do.fast = FALSE) 
cellChat <- identifyOverExpressedInteractions(cellChat)

# The number of highly variable ligand-receptor pairs used for signaling inference is <fill in here>

## 6. Calculate communication probability and infer cellular communication network
# In this step, CellChat models cell-cell communication through network analysis (creating a network on a per-pathway basis)
# In this network: nodes are celltypes and edges are the likelihood that different cell types are interacting (communication probability)
# From these networks, cell-cell communication probabilities can be inferred 

# Here we are setting the "type" to be based off the "triMean" for a given sample
# By setting the interaction type to "triMean" we are looking for fewer but overall stronger cell-cell communication interactions 
# "triMean" is not the only one way to run this, see the CellChat GitHub/vignette for more ways to run this command

# Note** This step takes a very long time to run on the AWS because the sever is pretty small
# Because of this we have a pre-made cellChat object in the cookingShow/ directory you will use for the next worksheet
# However - in real life you would run the command below (commented out because we are not running it today):

#cellChat <- computeCommunProb(cellChat, type = "triMean")





