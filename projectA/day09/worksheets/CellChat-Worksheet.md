# Getting started with CellChat
### Author: Georgia Barone (2024)

Today, we will explore one way to analyze processed scRNA-seq data with the package CellChat. CellChat is a networks-based R package that infers the likelihood of two cell types communicating. This package quantifies the signaling communication probability between two cell types by assessing the expression of ligand-receptor pairs within a pre-labeled scRNA-seq Seurat object.

## 1. Log on to RStudio 
If you successfully installed Seurat & CellChat locally, pull up your local RStudio. If you could not get CellChat installed locally, log onto the RStudio on the AWS server. 

If you are on the AWS server, please select R version 4.4.0

![AWS-viz-R-version.png](./day9-screenshots/AWS-viz-R-version.png)

## 2. Complete R scripts for today
Complete the scripts in the `/Users/<your-username>/srworkshop/projectA/day09/scripts/` directory.

Start with: `Cell-Chat-Part-1.R`

When you finish the first worksheet, complete: ``Cell-Chat-Part-2.R`

**NOTE: If you are running CellChat locally, you will need to transfer the following R objects from the AWS server to your local computer to complete the workheets and analysis for today**

Path to Seurat object for `Cell-Chat-Part-1.R`:
- Trisomy 21 Seurat object (what to start with):
`/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/t21_official_umap_clust_06.25.24.Rdata`
- Disomy 21 Seurat object (if you finish early, feel free to complete the worksheets a second time with the Disomy 21 objects):`/scratch/Shares/public/sread2024/cookingShow/day8a/labeled-seurat-objs/d21_subset_official_umap_clust_06.25.24.Rdata`

Path to CellChat object for `Cell-Chat-Part-2.R`:
- Trisomy 21 CellChat object (what to start with): `/scratch/Shares/public/sread2024/cookingShow/day9a/cellchat/cellchat_object-t21.RData`
- Disomy 21 CellChat object (for if you finish early): `/scratch/Shares/public/sread2024/cookingShow/day9a/cellchat/cellchat_object-d21.RData`




