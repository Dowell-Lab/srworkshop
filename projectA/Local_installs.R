## CRAN packages
install.packages(c(
  "Seurat",
  "patchwork",
  "dplyr",
  "ggplot2",
  "Matrix",
  "R.utils"
))

## Bioconductor packages (SingleR, celldex, monocle3 dependencies)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "SingleR",
  "celldex"
))

## monocle3 (from Bioconductor)
BiocManager::install("monocle3")

## SeuratWrappers (GitHub)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("satijalab/seurat-wrappers")

## CellChat (GitHub; requires devtools/remotes)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("sqjin/CellChat")


library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)
library(ggplot2)
library(Matrix)
library(monocle3)
library(R.utils)
library(SeuratWrappers)
