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
# -----------------------------------------------------------------------------
# HOW TO USE THIS TEMPLATE
#   Work through it top to bottom alongside the worksheet (05_pseudotime.md).
#   The boilerplate is already written for you. Wherever you see a block like:
#
#       # ---- Step 3a: preprocess ... ----
#       # Hint: preprocess_cds(); args ...
#       # YOUR CODE HERE:
#       cds <-
#
#   finish the line yourself. If you get stuck, the completed answer key is in
#   scripts_finished/05_pseudotime.R -- try it on your own first.
# =============================================================================

source("~/srworkshop/projectA/00_paths_and_setup.R")

library(monocle3)
library(R.utils)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
#If you didn't make this yet, use mine!
#combined <- readRDS(file.path(COOKING, "gata1_combined_annotated_joined.rds"))


# ---- 1. Inspect the lineage labels (already written) ------------------------
# A trajectory only makes sense within a single continuum. Erythroid maturation
# is the DOMINANT axis in this dataset, so we run monocle3 on all the cells here.
# Print the SingleR labels first so you can see what your run produced (and, if
# you wanted to, which pattern you would subset on).
print(sort(table(combined$SingleR_labels_other), decreasing = TRUE))


# ---- 2. Convert Seurat -> monocle3 cell_data_set ----------------------------
# ---- Step 2a: build the cell_data_set (cds) ---------------------------------
# monocle3 works on its own object type. Convert the Seurat object into a cds.
# Hint: SeuratWrappers::as.cell_data_set(); args combined
# YOUR CODE HERE:
cds <-
#Ignore the warning about cluser patriations. We will run cluster_cells in a minute

# Carry gene names over so plot_cells(genes = ...) can find them. (written for you)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(combined[["RNA"]])

# =============================================================================
# BASIC pseudotime - required setup + ordering
# =============================================================================
# Required order: preprocess -> reduce_dimension -> cluster -> learn_graph
# monocle3 is strict about this order -- run the four calls below in sequence.

# ---- Step 3a: preprocess the cds (PCA, 10 dims) -----------------------------
# PCA on the cds, matching the 10 PCs we used back in script 02.
# Hint: monocle3::preprocess_cds(); args cds, method = "PCA", num_dim = 10
# YOUR CODE HERE:
cds <-

plot = monocle3::plot_pc_variance_explained(cds)
# you can ignore this warning, monocole has not updated.
plot

ggsave(file.path(OUT_DIR, "gata1_pt_pc_variance.png"),
       plot, width = 6, height = 4, dpi = 150)

# ---- Step 3b: reduce dimensions (UMAP) --------------------------------------
# Embed the preprocessed cds into 2D so the trajectory can be drawn.
# Hint: monocle3::reduce_dimension(); args cds, reduction_method = "UMAP",
#       preprocess_method = "PCA"
# YOUR CODE HERE:
cds <-

# ---- Step 3c: cluster the cells ---------------------------------------------
# monocle3 finds its own clusters/partitions; learn_graph needs these.
# Hint: monocle3::cluster_cells(); args cds
# YOUR CODE HERE:
cds <-

# ---- Step 3d: learn the principal graph -------------------------------------
# Fit the principal graph -- the "backbone" the trajectory follows.
# Hint: monocle3::learn_graph(); args cds
# YOUR CODE HERE:
cds <-
#ignore the warning

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

# ---- Colour the graph by day (WORKED EXAMPLE) -------------------------------
# This first plot_cells is done for you. Notice the pattern:
#   plot_cells(cds, color_cells_by = <a metadata column>, label_cell_groups = FALSE)
# then print and save it. Copy this pattern for the one below.
plot <- monocle3::plot_cells(cds, color_cells_by = "day", label_cell_groups = FALSE)
plot
save_graph(plot,
           "gata1_pt_by_day.png")

# ---- Step 3e: colour the graph by construct ---------------------------------
# Same plot_cells pattern, coloured by "construct" instead of "day".
# Hint: monocle3::plot_cells(); args cds, color_cells_by = "construct",
#       label_cell_groups = FALSE
# YOUR CODE HERE:
plot <-
plot
save_graph(plot,
           "gata1_pt_by_construct.png")

# Erythroid maturation markers to orient the path (early -> late). (written for you)
plot <- monocle3::plot_cells(cds, genes = c("CD34", "GATA1", "TFRC", "HBG1"),label_cell_groups = FALSE, show_trajectory_graph = FALSE)
plot
save_graph(plot,
           "gata1_pt_markers.png", w = 9, h = 7)

# ---- 3. Pick the root WITHOUT a mouse click ---------------------------------
# order_cells() normally opens an interactive node-picker. For a headless run we
# pick the root automatically as the node sitting among the EARLIEST cells -
# here, the cells from day D7 (progenitor-rich, high CD34). The helper below is
# written for you.
get_earliest_node <- function(cds, time_col = "day", early_value = "D7") {
  cell_ids <- which(colData(cds)[[time_col]] == early_value)
  if (length(cell_ids) == 0) return(NULL)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  top_node <- as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[top_node]
}

# ---- Step 4a: find the earliest (D7) node -----------------------------------
# Use the helper to grab the graph node sitting among the D7 cells.
# Hint: get_earliest_node(); args cds, "day", "D7"
# YOUR CODE HERE:
root_node <-

if (!is.null(root_node)) {
  # ---- Step 4b: order the cells from that root ------------------------------
  # Set pseudotime = 0 at root_node and order every cell along the graph.
  # Hint: monocle3::order_cells(); args cds, root_pr_nodes = root_node
  # YOUR CODE HERE:
  cds <-
} else {
  # Fallback: interactive picker (RStudio only). (written for you)
  cds <- monocle3::order_cells(cds)
}

save_graph(monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                                label_cell_groups = FALSE, label_leaves = FALSE,
                                label_branch_points = FALSE, graph_label_size = 1.5),
           "gata1_pt_pseudotime.png")

# ---- Step 5a: extract the pseudotime values ---------------------------------
# VALIDATE: does pseudotime track the experimental day? It should, broadly.
# Pull the per-cell pseudotime out of the cds so we can plot it against day.
# Hint: monocle3::pseudotime(); args cds, reduction_method = "UMAP"
# YOUR CODE HERE:
pt <-

# Build the validation frame and violin plot. (written for you)
val <- data.frame(pseudotime = pt, day = colData(cds)$day)
val <- val[is.finite(val$pseudotime), ]
p_val <- ggplot(val, aes(x = day, y = pseudotime, fill = day)) +
  geom_violin() + geom_boxplot(width = 0.1, outlier.size = 0.3) +
  labs(title = "Pseudotime vs. experimental day (sanity check)") + theme_bw()
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

# ---- Step 6a: test which genes vary along the graph -------------------------
# Run monocle3's spatial test over the principal graph.
# Hint: monocle3::graph_test(); args cds, neighbor_graph = "principal_graph",
#       cores = 4
# YOUR CODE HERE:
cds_pr_test_res <-

# ---- Step 6b: keep the strong, significant trajectory genes -----------------
# Filter to genes with a clean gradient (high morans_I) that are significant.
# Hint: subset(); args cds_pr_test_res, q_value < 0.0005 & morans_I > 0.25
# YOUR CODE HERE:
deg <-
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
