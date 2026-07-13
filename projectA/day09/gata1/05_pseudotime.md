# GATA1 Lesson 05 — Trajectory analysis / pseudotime with monocle3

**Template script (you fill this in):** `05_pseudotime.R`
**Answer key:** `scripts_finished/05_pseudotime.R`
**Input:** `OUT_DIR/gata1_combined_annotated_joined.rds` (from Lesson 03)
**Output:** `OUT_DIR/gata1_cds_pseudotime.rds`

## How to use this worksheet

Open the template `05_pseudotime.R` in RStudio next to this worksheet. We work
through the script one section at a time. For each step you will see:

1. what the step does and why,
2. the code to write, and
3. a task cue indicating which `# YOUR CODE HERE:` block in the template to
   complete.

The boilerplate (loading libraries and data, the `save_graph` helper, the
`get_earliest_node` root-picker, the validation violin, and the advanced gene
plots) is already written for you. You complete the key line or lines in each
step. Attempt each one before consulting the answer key in `scripts_finished/`.

Each step provides the function name and its key arguments, so you assemble the
call rather than write it from scratch.

> Do not edit paths in the script. All paths live in
> `~/srworkshop/projectA/00_paths_and_setup.R`, and everything you create is
> written to `OUT_DIR`, not the read-only share.

## Overview

Lesson 03 gave every cell a SingleR label. Clusters and labels are *snapshots*;
development is a *continuum*. **Pseudotime** orders cells along an inferred
developmental path so you can watch a program unfold gene by gene
([Trapnell et al., *Nat Biotechnol* 2014](https://www.nature.com/articles/nbt.2859);
[Cao et al., *Nature* 2019](https://www.nature.com/articles/s41586-019-0969-x)).
The overall flow is:

> inspect labels → Seurat→cds → preprocess → reduce → cluster → learn_graph →
> pick root → order cells → validate vs. day → trajectory genes

> **Why this dataset is ideal for teaching pseudotime.** GSE271399 is a
> differentiation **time course** — cells were collected at D7, D9, and D11. So
> we have an independent, experimental clock to check the algorithm against.
> Pseudotime is *inferred* from expression alone; if it broadly agrees with the
> real sampling day, you can trust it. That validation step (Step 5) is the heart
> of this lesson.

---

## Step 0 — Setup and load (already written)

There is nothing to complete here; read and run this block. The script sources
the shared paths file, loads the libraries (note the new ones: `monocle3`,
`SeuratWrappers`, `R.utils`), and reads the annotated object saved by Lesson 03.

```r
source("~/srworkshop/projectA/00_paths_and_setup.R")
library(monocle3)
library(R.utils)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
```

> If you just finished Lesson 03 in the same R session, the annotated object is
> already in your environment. Otherwise reload it as above — or, if you never
> made it, use the pre-baked copy in `COOKING` (the commented `readRDS` line right
> below it in the script).

---

## Step 1 — Inspect the lineage labels (already written)

A trajectory only makes sense within a single continuum — you can't order
unrelated cell types along one axis. Erythroid maturation is the **dominant**
process in this dataset, so here we run monocle3 on all the cells and let the
trajectory follow that axis. Before doing anything, print the SingleR labels so
you can see what your run actually produced (and which pattern you *would* subset
on if a run were dominated by several lineages):

```r
print(sort(table(combined$SingleR_labels_other), decreasing = TRUE))
```

This block is written for you — just read the table. If a future dataset weren't
so cleanly erythroid, this is where you'd `subset()` to a single lineage first
(see the Common pitfalls).

---

## Step 2 — Seurat → monocle3

`SeuratWrappers::as.cell_data_set()` converts the Seurat object into a monocle3
`cell_data_set` (cds) — monocle3 works on its own object type, not on a Seurat
object.

```r
cds <- SeuratWrappers::as.cell_data_set(combined)
```

You can ignore the warning about cluster partitions — we run `cluster_cells()` in
Step 3. The next line (written for you) copies gene names into `gene_short_name`
so that `plot_cells(genes = ...)` can find them by symbol.

**Your task:** Complete **Step 2a** in the template.

---

## Step 3 — The required monocle3 order

monocle3 has a **fixed** pipeline; run these four calls in this order or it
errors / gives garbage:

1. `preprocess_cds(method = "PCA", num_dim = 10)` — PCA, matching the 10 PCs we
   chose in Lesson 02.
2. `reduce_dimension(reduction_method = "UMAP")` — the 2D embedding.
3. `cluster_cells()` — monocle3's own partitions (`learn_graph` needs them).
4. `learn_graph()` — fits the principal graph (the "backbone" the trajectory
   follows).

```r
cds <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = 10)
cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds)
```

**Your task:** Complete **Steps 3a, 3b, 3c, and 3d** in the template — one call
each, in order. The PC-variance plot between 3a and 3b is written for you.

Then we plot the graph colored by what we know about each cell. The `day` plot is
provided as a **worked example**; you repeat the same `plot_cells` pattern for
`construct`, and the erythroid-marker plot (`CD34` early → `GATA1`/`TFRC` →
`HBG1` late, to orient early vs. late) is written for you.

```r
plot <- monocle3::plot_cells(cds, color_cells_by = "construct", label_cell_groups = FALSE)
```

**Your task:** Complete **Step 3e** in the template.

---

## Step 4 — Choosing the root without a mouse click

`order_cells()` normally opens an **interactive** window to click the starting
node — fine in RStudio, impossible on a headless server. So the script includes
`get_earliest_node()` (written for you), which picks the graph node sitting among
the earliest cells — the D7, progenitor-rich cells — automatically. If that helper
can't find D7 cells it falls back to the interactive picker. Either way, the root
defines pseudotime = 0.

```r
root_node <- get_earliest_node(cds, "day", "D7")
if (!is.null(root_node)) {
  cds <- monocle3::order_cells(cds, root_pr_nodes = root_node)
} else {
  cds <- monocle3::order_cells(cds)   # interactive fallback (written for you)
}
```

**Your task:** Complete **Step 4a** (find the node) and **Step 4b** (order the
cells from it) in the template.

---

## Step 5 — Validate pseudotime against day

This is the payoff. We extract `pseudotime(cds)` and plot it against the
experimental `day`:

```r
pt <- monocle3::pseudotime(cds, reduction_method = "UMAP")
```

The data frame and violin plot that follow are written for you:

```r
ggplot(val, aes(x = day, y = pseudotime, fill = day)) + geom_violin() ...
```

**Your task:** Complete **Step 5a** in the template (extract `pt`).

If D7 cells sit at low pseudotime and D11 cells at high pseudotime, the inferred
ordering recovered the real biology — strong evidence the trajectory is
trustworthy. If they're scrambled, something is wrong (wrong root, or the cells
aren't really one continuum). **Always validate pseudotime against something
external when you can.** The script then saves the cds to
`gata1_cds_pseudotime.rds` (written for you).

---

## Step 6 — Genes that change along the trajectory (advanced)

`graph_test(neighbor_graph = "principal_graph")` finds genes whose expression
varies *spatially* along the trajectory:

- **`morans_I`** — spatial autocorrelation; higher means a cleaner gradient
  along the path.
- **`q_value`** — multiple-testing-corrected significance.

```r
cds_pr_test_res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
deg <- subset(cds_pr_test_res, q_value < 0.0005 & morans_I > 0.25)
```

We keep the strong, significant genes (`q_value < 0.0005 & morans_I > 0.25`), save
the gene list, and plot key erythroid genes with `plot_genes_in_pseudotime()` to
watch the maturation program turn on and off in order (that plot is written for
you). `choose_graph_segments()` (interactive) can isolate a single branch for
finer analysis.

**Your task:** Complete **Steps 6a and 6b** in the template.

---

## Common pitfalls

- **Running pseudotime on a mix of unrelated cell types.** The trajectory becomes
  meaningless. This dataset is erythroid-dominated so we run on all cells, but in
  a mixed dataset you'd `subset()` to one lineage first (Step 1 is where you'd
  check).
- **Picking the wrong root.** Pseudotime direction is arbitrary until you set the
  root; anchor it with biology (progenitors / earliest day / `CD34` high).
- **Not validating.** Pseudotime always produces *an* ordering. Check it against
  day, known markers, or a held-out signal before believing it.
- **Skipping steps.** monocle3's `preprocess → reduce → cluster → learn_graph`
  order is mandatory; calling them out of order errors or gives garbage.
- **Interactive calls on a server.** `order_cells()`/`choose_graph_segments()`
  need a display. Use the programmatic root helper for headless runs.

## Check your understanding

1. Why does a trajectory only make sense within a single continuum, and how would
   you subset to one lineage if a dataset had several?
2. How does the experimental `day` metadata let you validate pseudotime, and
   what would a *failed* validation look like?
3. What determines the direction of pseudotime, and how does the script set it
   without a mouse click?
4. What do `morans_I` and `q_value` measure in `graph_test`, and why use both?
5. Your pseudotime violin shows D11 cells at *lower* pseudotime than D7. List two
   things you'd check.
