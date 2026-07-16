# =============================================================================
# GATA1 track - Script 05: Trajectory / pseudotime with monocle3
# -----------------------------------------------------------------------------
# Orders cells along a developmental continuum. This dataset is a hematopoietic
# differentiation TIME COURSE (collected at D7, D9, D11), so pseudotime has a
# real-world anchor: the inferred ordering should broadly track the sampling day.
# That makes it a great teaching case - you can VALIDATE pseudotime against the
# experimental day metadata.
#
# Input : OUT_DIR/gata1_combined_annotated_joined.rds  (from script 03)
# Output: a monocle3 cds + trajectory plots in OUT_DIR.
#
# Interactive calls (order_cells / choose_cells) pop up a chooser window when
# run in RStudio. A non-interactive root-picking helper is provided so the
# script also runs start-to-finish on a headless server.
#
# ANSWER KEY for the fill-in template ../05_pseudotime.R. Try the template
# (with the worksheet 05_pseudotime.md) before reading this.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(monocle3)
library(R.utils)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

#combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined_subsampled.rds"))
#If you didn't make this yet, use mine!
combined <- readRDS(file.path(COOKING, "gata1_combined_annotated_joined_subsampled.rds"))


# ---- 1. Inspect the lineage labels (already written) ------------------------
# A trajectory only makes sense within a single continuum. Erythroid maturation
# is the DOMINANT axis in this dataset, so we run monocle3 on all the cells here.
# Print the SingleR labels first so you can see what your run produced (and, if
# you wanted to, which pattern you would subset on).
#base R version
print(sort(table(combined$SingleR_labels_other), decreasing = TRUE))
#dplyr version
#as.data.frame(combined@meta.data) %>% group_by(SingleR_labels_other) %>% tally() %>% arrange(desc(n))


# ---- 2. Convert Seurat -> monocle3 cell_data_set ----------------------------
cds <- SeuratWrappers::as.cell_data_set(combined)
#Ignore the warning about cluser patriations. We will run cluster_cells in a minute

# Carry gene names over so plot_cells(genes = ...) can find them.
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(combined[["RNA"]])

# =============================================================================
# BASIC pseudotime - required setup + ordering
# =============================================================================
# Required order: preprocess -> reduce_dimension -> cluster -> learn_graph
cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = 6)

plot = monocle3::plot_pc_variance_explained(cds)
# you can ignore this warning, monocole has not updated. 
plot

ggsave(file.path(OUT_DIR, "gata1_pt_pc_variance.png"),
       plot, width = 6, height = 4, dpi = 150)

cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds)
#ignore the warning


#alternative to 2 keeping the umap and pca from seurat
#in the code with seurat (02 lesson)
#seurat_umap <- Embeddings(seurat_obj, reduction = "umap")    # matrix, cells × dims [web:51]
#seurat_pca <- Embeddings(seurat_obj, reduction = "pca")
#
#saveRDS(seurat_umap, file = "seurat_umap_coords.rds")
#saveRDS(seurat_pca,  file = "seurat_pca_coords.rds")

# Later, in monocole run

# Create cds_new from counts/metadata, then inject:
#seurat_umap  <- readRDS("seurat_umap_coords.rds")
#seurat_pca   <- readRDS("seurat_pca_coords.rds")
#combined <- readRDS(file.path(COOKING, "gata1_combined_annotated_joined_subsampled.rds"))
#cds_new <- SeuratWrappers::as.cell_data_set(combined)
#cds_new@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(combined[["RNA"]])
#reducedDims(cds_new)$UMAP <- seurat_umap
#reducedDims(cds_new)$PCA  <- seurat_pca


p<- plot_cells(
  cds,
  color_cells_by = "cluster",   # or "partition", "cell_type", etc.
  label_groups_by_cluster = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

p


# Color the trajectory by the things we know about each cell.
save_graph <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, file), plot = p, width = w, height = h, dpi = 150)
}

plot <- monocle3::plot_cells(cds, color_cells_by = "day", label_cell_groups = FALSE)
plot
save_graph(plot,
           "gata1_pt_by_day.png")
plot <- monocle3::plot_cells(cds, color_cells_by = "construct", label_cell_groups = FALSE)
plot
save_graph(plot,
           "gata1_pt_by_construct.png")

plot <- monocle3::plot_cells(cds, genes = c("CD34", "GATA1", "TFRC", "HBG1"),label_cell_groups = FALSE, show_trajectory_graph = FALSE)
plot
# Erythroid maturation markers to orient the path (early -> late).
save_graph(plot,
           "gata1_pt_markers.png", w = 9, h = 7)

# ---- 3. Pick the root WITHOUT a mouse click ---------------------------------
# order_cells() normally opens an interactive node-picker. For a headless run we
# pick the root automatically as the node sitting among the EARLIEST cells -
# here, the cells from day D7 (progenitor-rich, high CD34).
get_earliest_node <- function(cds, time_col = "day", early_value = "D7") {
  cell_ids <- which(colData(cds)[[time_col]] == early_value)
  if (length(cell_ids) == 0) return(NULL)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  top_node <- as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[top_node]
}

root_node <- get_earliest_node(cds, "day", "D7")
if (!is.null(root_node)) {
  cds <- monocle3::order_cells(cds, root_pr_nodes = root_node)
} else {
  # Fallback: interactive picker (RStudio only).
  cds <- monocle3::order_cells(cds)
}



p<- monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                         label_cell_groups = FALSE, label_leaves = FALSE,
                         label_branch_points = FALSE, graph_label_size = 1.5)
p

save_graph(p,
           "gata1_pt_pseudotime.png")

#cds <- monocle3::order_cells(cds)
cds_sub1 <- choose_graph_segments(cds, clear_cds = FALSE)
cds_sub2 <- choose_graph_segments(cds, clear_cds = FALSE)

p<- monocle3::plot_cells(cds_sub1, color_cells_by = "pseudotime",
                         label_cell_groups = FALSE, label_leaves = FALSE,
                         label_branch_points = FALSE, graph_label_size = 1.5)
p

save_graph(p,
           "gata1_pt_pseudotime_branch1.png")

p<- monocle3::plot_cells(cds_sub1, color_cells_by = "pseudotime")

# VALIDATE: does pseudotime track the experimental day? It should, broadly.
pt <- monocle3::pseudotime(cds, reduction_method = "UMAP")
val <- data.frame(pseudotime = pt, day = colData(cds)$day)
val <- val[is.finite(val$pseudotime), ]
p_val <- ggplot(val, aes(x = day, y = pseudotime, fill = day)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.size = 0.3) +
  labs(title = "Pseudotime vs. experimental day (sanity check)") + theme_bw()
p_val
save_graph(p_val, "gata1_pt_vs_day_violin.png")


saveRDS(cds, file.path(OUT_DIR, "gata1_cds_pseudotime.rds"))
message("Saved: ", file.path(OUT_DIR, "gata1_cds_pseudotime.rds"))



#cds<- readRDS(file.path(OUT_DIR, "gata1_cds_pseudotime.rds"))
# Ignore the warning about save_monocle_objects()
#it quite working when they update one part and save_monocle_objects wont work with it.

# =============================================================================
# ADVANCED - genes that change along the trajectory
# =============================================================================
# graph_test: which genes vary SPATIALLY along the trajectory graph?
#   morans_I : strength of spatial autocorrelation (higher = cleaner gradient).
#   q_value  : multiple-testing-corrected significance.
cds_pr_test_res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
deg <- subset(cds_pr_test_res, q_value < 0.0005 & morans_I > 0.25)
deg <- deg[order(-deg$morans_I), ]
pr_deg_ids <- rownames(deg)
print(head(pr_deg_ids, 30))
write.csv(deg, file.path(OUT_DIR, "gata1_pt_trajectory_genes.csv"))

# Plot a few erythroid genes ALONG pseudotime to see the maturation program.
genes_look_interesting <- intersect(c("CD34", "GATA1", "TFRC", "HBG1", "HBG2", "GFI1B"),
                                     rownames(rowData(cds)))
if (length(genes_look_interesting) > 0) {
  subset_cds <- cds[rowData(cds)$gene_short_name %in% genes_look_interesting, ]
  p_genes <- monocle3::plot_genes_in_pseudotime(subset_cds, color_cells_by = "day",
                                                vertical_jitter = TRUE, cell_size = 0.3)
  save_graph(p_genes, "gata1_pt_genes_in_pseudotime.png", w = 7, h = 8)
}

# choose_graph_segments() (interactive) isolates one branch start->end:
# cds_sub <- monocle3::choose_graph_segments(cds)
