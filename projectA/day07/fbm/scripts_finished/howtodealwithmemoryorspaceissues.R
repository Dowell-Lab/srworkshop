

# Option 1: subsample (only really works for class)
meta <- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")   # cell names become a column

cells_to_keep <- meta %>%
  group_by(day_construct_genotype) %>%
  slice_sample(n = 2000, replace = FALSE) %>%
  pull(cell_id)

# Subset Seurat object
combined_subsampled <- subset(obj, cells = cells_to_keep)


#option 2 parralelize or don't parralelelize 

library(future)
options(future.globals.maxSize = 12 * 1024^4)
#faster but takes more memory
plan(multisession, workers = 4)
#slower but takes less memory
plan(sequential)


# option 3: DietSeurat to keep only essential parts 

# 1. Drop scale.data once PCA/UMAP/clustering are done 

obj@assays$RNA@scale.data <- matrix() 
  if ("SCT" %in% names(merged@assays)) { merged@assays$SCT@scale.data <- matrix() } 

obj_slim <- DietSeurat( object = obj, assays = c("RNA", "SCT"), dimreducs = c("pca", "umap"), graphs = NULL, misc = TRUE ) 


# Also, often run rm the old object and 
gc() 
rm(merged) 
gc()