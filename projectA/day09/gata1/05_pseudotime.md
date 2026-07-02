# GATA1 Lesson 05 — Trajectory analysis / pseudotime with monocle3

**Script:** `scripts/gata1/05_pseudotime.R`
**Input:** `OUT_DIR/gata1_combined_annotated.rds` (from Lesson 03)

Clusters are snapshots; development is a continuum. **Pseudotime** orders cells
along an inferred developmental path so you can watch a program unfold gene by
gene ([Trapnell et al., *Nat Biotechnol* 2014](https://www.nature.com/articles/nbt.2859);
[Cao et al., *Nature* 2019](https://www.nature.com/articles/s41586-019-0969-x)).

> **Why this dataset is ideal for teaching pseudotime.** GSE271399 is a
> differentiation **time course** — cells were collected at D7, D9, and D11. So
> we have an independent, experimental clock to check the algorithm against.
> Pseudotime is *inferred* from expression alone; if it broadly agrees with the
> real sampling day, you can trust it. That validation step is the heart of this
> lesson.

## Step 1 — Subset to one lineage

A trajectory only makes sense within a single continuum — you can't order
unrelated cell types along one axis. Erythroid maturation is the dominant
process here, so we keep the erythroid-labeled cells from SingleR (matching
`"Erythro|Eryth|Ery"`). The script prints the label table first so you can
adjust the pattern to whatever your run produced.

## Step 2 — Seurat → monocle3

`SeuratWrappers::as.cell_data_set()` converts the Seurat object into a monocle3
`cell_data_set` (cds). We copy gene names into `gene_short_name` so that
`plot_cells(genes = ...)` can find them by symbol.

## Step 3 — The required monocle3 order

monocle3 has a fixed pipeline; run it in this order:

1. `preprocess_cds(method = "PCA", num_dim = 10)` — PCA, matching our 10 PCs.
2. `reduce_dimension(reduction_method = "UMAP")` — 2D embedding.
3. `cluster_cells()` — monocle3's own partitions.
4. `learn_graph()` — fits the principal graph (the "backbone" the trajectory
   follows).

We then plot the graph colored by `day`, by `construct`, and by erythroid
markers (`CD34` early → `GATA1`/`TFRC` → `HBG1` late) to orient early vs. late.

## Step 4 — Choosing the root without a mouse click

`order_cells()` normally opens an **interactive** window to click the starting
node — fine in RStudio, impossible on a headless server. So the script includes
`get_earliest_node()`, which picks the graph node sitting among the earliest
cells (the D7, progenitor-rich cells) automatically. If that helper can't find
D7 cells it falls back to the interactive picker. Either way, the root defines
pseudotime = 0.

## Step 5 — Validate pseudotime against day

This is the payoff. We extract `pseudotime(cds)` and plot it against the
experimental `day`:

```r
ggplot(val, aes(x = day, y = pseudotime, fill = day)) + geom_violin() ...
```

If D7 cells sit at low pseudotime and D11 cells at high pseudotime, the inferred
ordering recovered the real biology — strong evidence the trajectory is
trustworthy. If they're scrambled, something is wrong (wrong lineage subset,
wrong root, or not a real continuum). **Always validate pseudotime against
something external when you can.**

## Step 6 — Genes that change along the trajectory (advanced)

`graph_test(neighbor_graph = "principal_graph")` finds genes whose expression
varies *spatially* along the trajectory:

- **`morans_I`** — spatial autocorrelation; higher means a cleaner gradient
  along the path.
- **`q_value`** — multiple-testing-corrected significance.

We keep `q_value < 0.0005 & morans_I > 0.25`, save the gene list, and plot key
erythroid genes with `plot_genes_in_pseudotime()` to see the maturation program
turn on and off in order. `choose_graph_segments()` (interactive) can isolate a
single branch for finer analysis.

## Common pitfalls

- **Running pseudotime on all cell types at once.** The trajectory becomes
  meaningless. Subset to one lineage first.
- **Picking the wrong root.** Pseudotime direction is arbitrary until you set the
  root; anchor it with biology (progenitors / earliest day / `CD34` high).
- **Not validating.** Pseudotime always produces *an* ordering. Check it against
  day, known markers, or a held-out signal before believing it.
- **Skipping steps.** monocle3's `preprocess → reduce → cluster → learn_graph`
  order is mandatory; calling them out of order errors or gives garbage.
- **Interactive calls on a server.** `order_cells()`/`choose_graph_segments()`
  need a display. Use the programmatic root helper for headless runs.

## Check your understanding

1. Why must you subset to a single lineage before fitting a trajectory?
2. How does the experimental `day` metadata let you validate pseudotime, and
   what would a *failed* validation look like?
3. What determines the direction of pseudotime, and how does the script set it
   without a mouse click?
4. What do `morans_I` and `q_value` measure in `graph_test`, and why use both?
5. Your pseudotime violin shows D11 cells at *lower* pseudotime than D7. List two
   things you'd check.
