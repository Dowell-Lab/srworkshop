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

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))
install.packages("devtools")


remotes::install_github("bnprks/BPCells/r")
#alternativly

#If that doesn't work, you might have to install hdf5 on your computer, then go back and install the above. 
#CONDA_SUBDIR=osx-arm64 conda create -n osx-arm hdf5
#conda activate osx-arm
#/usr/local/bin/R

devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('immunogenomics/presto')
#alternativly
#pak::pak('cole-trapnell-lab/monocle3')


#another option if devtools does not work is pak::pak('your_package')


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
#alternativly
#pak::pak('sqjin/CellChat')


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
