# Lesson 5 — Pseudotime / Trajectory Analysis (monocle3)

**Script:** `scripts/05_pseudotime.R`
**Dataset:** the pre-labeled **T21** (Trisomy 21) fetal bone marrow object, subset to the **erythroid** lineage.
**Goal:** Order cells along a developmental continuum and find the genes that change as cells mature.

A UMAP is a snapshot. But many biological processes — differentiation, maturation — are *continuous*. **Pseudotime** assigns each cell a position along an inferred trajectory, so you can study a process you only ever sampled as a frozen mixture. We use erythroid (red blood cell) maturation because it's a clean, well-understood progression.

---

## Setup

### 0–1. Load and subset to one lineage
We load the labeled object and **subset to `erythroid`** using `Idents()` + `subset()`. Two reasons: a trajectory only makes sense within a lineage (you can't order an erythroid cell relative to a T cell), and a smaller object keeps the computation fast. `cell.labels` then gives the fine erythroid sub-stages (early/mid/late), which are the "ground truth" we hope pseudotime recovers.

### 2. Convert Seurat → monocle3 (`as.cell_data_set`)
monocle3 works on its own **`cell_data_set` (cds)** object, not a Seurat object. `SeuratWrappers::as.cell_data_set()` converts it. We then manually copy gene names into `gene_short_name` so that `plot_cells(genes = ...)` can find genes by symbol — a classic gotcha if you skip it.

---

## BASIC tier — the minimum to get a pseudotime

monocle3 requires **four setup steps, in order**, before you can order cells:

1. **`preprocess_cds(method="PCA", num_dim=3)`** — PCA dimensionality reduction. `num_dim` is small here because the erythroid subset is simple; `plot_pc_variance_explained` tells you how many PCs you actually need.
2. **`reduce_dimension(reduction_method="UMAP")`** — monocle3's own UMAP (independent of Seurat's).
3. **`cluster_cells`** — groups cells and, importantly, identifies **partitions** (disconnected pieces of the data).
4. **`learn_graph`** — fits a **principal graph**: the skeletal "tree" through the data that the trajectory follows.

Then:
- **Color the graph** by `cell.labels`, `orig.ident`, genes, etc. (`plot_cells`). Coloring by maturation markers (`RUNX1`, `HBB`, `PCNA`) reveals *where* the trajectory starts — you need this before ordering.
- **`order_cells()`** — **interactive**: a window opens and you click the **root** (the earliest, least mature cells). Pseudotime is measured as graph distance from that root.
- **`plot_cells(color_cells_by="pseudotime")`** — the payoff: cells colored by how far along the trajectory they are.
- **`pseudotime(cds)`** — extracts the numeric values for downstream stats.

> **Conceptual point:** pseudotime is *relative ordering*, not real clock time. It tells you sequence and relative spacing, not minutes or days. The root choice defines direction — choose it from biology, not convenience.

---

## ADVANCED tier — what changes along the trajectory

### `use_partition = FALSE`
Re-running `learn_graph(use_partition = FALSE)` forces monocle3 to treat all erythroid cells as **one connected progression** rather than separate partitions — appropriate when you believe it's a single continuous lineage.

### `graph_test` — trajectory-variable genes
This asks: *which genes vary spatially along the principal graph?*
- **Moran's I** measures spatial autocorrelation — high means the gene's expression is organized along the trajectory rather than scattered randomly.
- **`q_value`** is the multiple-testing-corrected significance.
We keep genes with `q_value < 0.0005 & morans_I > 0.5` — strong, significant trajectory genes.

### Handling ribosomal genes
In erythroid maturation, ribosomal genes (`RPL*`, `RPS*`) dominate the hit list because protein synthesis shuts down as cells enucleate. They're real but often not what you care about, so we split them out with `grepl` to focus on the more interesting remainder — a practical data-wrangling move worth teaching.

### Visualizing dynamics
- **`plot_genes_by_group`** — dot plot of trajectory genes across erythroid sub-stages (`maximal_on_diag` orders them so the pattern reads cleanly down the diagonal).
- **`plot_genes_in_pseudotime`** — expression of selected genes *as a function of pseudotime*, for chosen sub-stages. This is the clearest view of a gene turning on/off during maturation.

### Picking a sub-trajectory
- **`choose_graph_segments()`** (commented, interactive) — click a start and end node to isolate one branch.
- **`fit_models(model_formula_str = "~cell.labels")`** (commented) — a regression alternative to `graph_test` for modeling expression vs. a covariate.

---

## Common pitfalls
- Running pseudotime on a mix of unrelated cell types — order is meaningless across lineages.
- Forgetting to set `gene_short_name`, then `plot_cells(genes=...)` can't find genes.
- Choosing the root by convenience instead of biology — it sets the entire direction.
- Over-interpreting pseudotime as real elapsed time.

## Check your understanding
1. Why must you subset to a single lineage before building a trajectory?
2. What does a high Moran's I tell you about a gene?
3. Why do ribosomal genes dominate the erythroid trajectory gene list, and why might you set them aside?
4. What does the root you click in `order_cells()` actually determine?
