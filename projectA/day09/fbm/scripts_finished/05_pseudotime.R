# =============================================================================
# 05_pseudotime.R
# -----------------------------------------------------------------------------
# Step 5: trajectory / pseudotime analysis with monocle3.
# We order cells along a developmental continuum (here: erythroid maturation).
#
# DATASET: the pre-labeled T21 (Trisomy 21) fetal bone marrow object. We subset
# to the erythroid lineage so the trajectory is biologically meaningful and the
# computation stays fast.
#
# This script has THREE tiers:
#   BASIC    : minimal required steps to get a pseudotime.
#   ADVANCED : trajectory-variable genes, segment selection, gene dynamics.
# (Interactive calls like order_cells()/choose_cells() pop up a chooser window.)
#
# Companion lesson: lessons/05_pseudotime.md
# =============================================================================

source("../00_paths_and_setup.R")   # gives T21_LABELED, load_one()

library(monocle3)
library(R.utils)
library(Seurat)
library(SeuratWrappers)

# -----------------------------------------------------------------------------
# 0. Load the labeled object and inspect it
# -----------------------------------------------------------------------------
t21 <- load_one(T21_LABELED)
ls()

# What cell types and samples are present?
unique(t21@meta.data$broad_extfig7A_cell.labels)
unique(t21@meta.data$orig.ident)
table(t21@meta.data$cell.labels)

# -----------------------------------------------------------------------------
# 1. Subset to ONE lineage (erythroid) so the trajectory is interpretable
# -----------------------------------------------------------------------------
# Set the active identity to the broad label column, then keep only "erythroid".
Idents(object = t21) <- "broad_extfig7A_cell.labels"
red_T21 <- subset(t21, idents = c("erythroid"), invert = FALSE)
table(red_T21@meta.data$cell.labels)   # the fine erythroid sub-stages

# -----------------------------------------------------------------------------
# 2. Convert the Seurat object to a monocle3 cell_data_set (cds)
# -----------------------------------------------------------------------------
cds <- SeuratWrappers::as.cell_data_set(red_T21)
# Carry over gene names so plot_cells(genes=...) can find them.
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(red_T21[["RNA"]])

# =============================================================================
# BASIC pseudotime — the four required setup steps + ordering
# =============================================================================
# Required before pseudotime: preprocess -> reduce_dimension -> cluster -> learn_graph
cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = 3)
monocle3::plot_pc_variance_explained(cds)   # how much variance each PC explains

cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds)           # fit the principal-graph trajectory

# What does it look like, and what can we color by?
monocle3::plot_cells(cds)
colnames(monocle3::pData(cds))
monocle3::plot_cells(cds, color_cells_by = "cell.labels")
monocle3::plot_cells(cds, color_cells_by = "orig.ident", label_cell_groups = FALSE)

# Color by genes that mark erythroid maturation to find the START of the path.
rownames(cds)[grepl("RUNX", rownames(cds))]
monocle3::plot_cells(cds, genes = c("RUNX1", "HBB", "PCNA"),
                     label_cell_groups = FALSE, show_trajectory_graph = FALSE)

# order_cells() opens an interactive window: click the root (earliest) node.
cds <- monocle3::order_cells(cds)
monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                     label_cell_groups = FALSE, label_leaves = FALSE,
                     label_branch_points = FALSE, graph_label_size = 1.5)

# Pull the numeric pseudotime values out.
pt <- monocle3::pseudotime(cds, reduction_method = "UMAP")

# =============================================================================
# ADVANCED — genes that change along the trajectory + segment selection
# =============================================================================
# Treat all erythroid cells as one continuous progression (use_partition=FALSE).
cds <- monocle3::learn_graph(cds, use_partition = FALSE)
cds <- monocle3::order_cells(cds, reduction_method = "UMAP")

# graph_test: which genes vary SPATIALLY along the trajectory graph?
#   Moran's I  : how strongly a gene's expression is spatially autocorrelated.
#   q_value    : multiple-testing-corrected significance.
cds_pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
deg <- subset(cds_pr_test_res, q_value < 0.0005 & morans_I > 0.5)
deg <- deg[order(deg$q_value), ]
pr_deg_ids <- rownames(deg)

# Many top hits are ribosomal (RPL*/RPS*); separate them from the rest.
ribo_L <- pr_deg_ids[grepl("^RPL", pr_deg_ids)]
ribo_S <- pr_deg_ids[grepl("^RPS", pr_deg_ids)]
pr_deg_ids_other <- setdiff(pr_deg_ids, c(ribo_L, ribo_S))

# Dot plot of the non-ribosomal trajectory genes grouped by erythroid sub-stage.
plot_genes_by_group(cds, pr_deg_ids_other, group_cells_by = "cell.labels",
                    ordering_type = "maximal_on_diag")

# Plot a handful of interesting genes ALONG pseudotime for two sub-stages.
genes_look_interesting <- c("HMGN2", "RPL22", "HBB", "BLVRB")
subset_cds <- cds[rowData(cds)$gene_short_name %in% genes_look_interesting,
                  colData(cds)$cell.labels %in% c("early erythroid", "mid erythroid")]
monocle3::plot_genes_in_pseudotime(subset_cds, vertical_jitter = TRUE, cell_size = 0.1)

# choose_graph_segments() lets you interactively pick a START and END node to
# isolate one branch of the trajectory.
# cds_sub <- choose_graph_segments(cds)
# monocle3::plot_cells(cds_sub, color_cells_by = "pseudotime",
#                      label_cell_groups = FALSE, label_leaves = FALSE,
#                      label_branch_points = FALSE, graph_label_size = 1.5)

# Regression-style modeling of expression vs. cell label (alternative to graph_test).
# gene_fits <- fit_models(cds, model_formula_str = "~cell.labels")
# fit_coefs <- coefficient_table(gene_fits)
