###
#MAKING A TINY OBJECT FOR CLASS
###

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(Seurat)
library(tidyverse)
library(future)
plan(sequential)


combined <- readRDS(file.path(OUT_DIR, "gata1_combined_clustered.rds"))

meta <- combined@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")   # cell names become a column

cells_to_keep <- meta %>%
  group_by(day_construct_genotype) %>%
  slice_sample(n = 1000, replace = FALSE) %>%
  pull(cell_id)

# Subset Seurat object
combined_subsampled <- subset(combined, cells = cells_to_keep)

saveRDS(combined_subsampled, file.path(OUT_DIR, "gata1_combined_clustered_subsampled.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_clustered_subsampled.rds"))

rm(combined)
gc()
combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))

meta <- combined@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")   # cell names become a column

cells_to_keep <- meta %>%
  group_by(day_construct_genotype) %>%
  slice_sample(n = 1000, replace = FALSE) %>%
  pull(cell_id)

# Subset Seurat object
combined_subsampled <- subset(combined, cells = cells_to_keep)

saveRDS(combined_subsampled, file.path(OUT_DIR, "gata1_combined_annotated_joined_subsampled.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_combined_annotated_joined_subsampled.rds"))
