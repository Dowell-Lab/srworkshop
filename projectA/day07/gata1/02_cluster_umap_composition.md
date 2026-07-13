# GATA1 Lesson 02 — Normalization, PCA, clustering, UMAP, and composition

**Template script (you fill this in):** `02_cluster_umap_composition.R`
**Answer key:** `scripts_finished/02_cluster_umap_composition.R`
**Input:** `OUT_DIR/gata1_combined_qc.rds` (from Lesson 01)

## How to use this worksheet

Open the template `02_cluster_umap_composition.R` in RStudio next to this
worksheet. We work through the script one section at a time. For each step you
will see:

1. what the step does and why,
2. the code to write, and
3. a task cue indicating which `# YOUR CODE HERE:` block in the template to
   complete.

The boilerplate (loading data, the `ggplot` heatmaps, saving plots) is already
written for you. You complete the key line or lines in each step. Attempt each
one before consulting the answer key in `scripts_finished/`.

Each step provides the function name and its key arguments, so you assemble the
call rather than write it from scratch.

## Overview

Lesson 01 produced a filtered object with clean metadata. In this lesson we turn
raw counts into a map of cell states (clusters on a UMAP) and then quantify how
the experimental conditions are distributed across that map. The overall flow is:

> normalize → variable genes → (cell cycle) → scale → PCA → UMAP + clusters →
> visualize → composition

---

## Step 0 — Setup and load (already written)

There is nothing to complete here; read and run this block. The script sources
the shared paths file, loads the libraries, reads the object saved by Lesson 01,
and sets the default assay to RNA.

```r
source("~/srworkshop/projectA/00_paths_and_setup.R")
library(Seurat)
library(dplyr)
library(ggplot2)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_qc.rds"))
DefaultAssay(combined) <- "RNA"
```

> If you just finished Lesson 01 in the same R session, `combined` is already in
> your environment — you don't need to reload it, so comment out the `readRDS`
> line. Only run it on a fresh session (or after clearing your workspace).

> Do not edit paths in the script. All paths live in `00_paths_and_setup.R`, and
> everything you create is written to `OUT_DIR`, not the read-only share.

---

## Step 1 — Normalize and select variable genes

**Normalize.** Cells differ in total counts largely for technical reasons, such
as how many transcripts were captured. `NormalizeData` (LogNormalize) divides
each count by the cell's total, scales to 10,000, and log-transforms the result
so that expression is comparable across cells.

```r
combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = TRUE)
```

**Your task:** Complete **Step 1a** in the template.

**Variable genes.** Most genes are uninformative housekeeping genes. We keep the
3,000 most variable genes, which carry most of the biological signal, using
`FindVariableFeatures`.

```r
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                 nfeatures = 3000, verbose = TRUE)
head(VariableFeatures(combined), 20)
```

**Your task:** Complete **Step 1b** in the template. The `head(...)` line that
prints the top 20 variable genes is already written for you.

---

## Step 2 — Cell-cycle scoring

Dividing cells can form clusters that reflect the cell cycle rather than the
biology of interest. Before scaling, we score each cell's position in the cell
cycle so that this effect is visible and can be removed later if needed.

`CellCycleScoring` takes two curated gene lists — S-phase and G2M-phase markers,
which Seurat provides in `cc.genes.updated.2019` — and adds three columns to the
metadata: `S.Score`, `G2M.Score`, and `Phase` (G1, S, or G2M).

The gene lists are provided:

```r
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
```

Complete the scoring call:

```r
combined <- CellCycleScoring(
  combined,
  s.features   = s.genes,
  g2m.features = g2m.genes,
  set.ident    = TRUE
)
```

**Your task:** Complete **Step 1c** in the template.

---

## Step 3 — Scale and regress out technical drivers

`ScaleData` centers and scales each gene to mean 0 and variance 1. It can also
regress out unwanted signals. We regress out sequencing depth (`nCount_RNA`) and
mitochondrial fraction (`percent.mt`) so that these technical drivers are not
mistaken for biological structure in the PCA.

```r
combined <- ScaleData(
  combined,
  features        = VariableFeatures(combined),
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose         = TRUE
)
```

> Heads up: with the regression, this is the **slow step** — it can take a while
> (up to ~20 minutes). Keep `verbose = TRUE` so you can watch its progress
> rather than staring at a frozen prompt.

**Your task:** Complete **Step 1d** in the template.

> Optional experiment: the template includes a commented-out `ScaleData` that
> also regresses out `S.Score` and `G2M.Score`. If your `Phase` UMAP (Step 6)
> shows clusters splitting by cell cycle, run that version and compare the
> results.

---

## Step 4 — PCA and choosing the number of PCs

`RunPCA` compresses 3,000 genes into a few dozen principal components.

```r
combined <- RunPCA(combined, features = VariableFeatures(combined),
                   npcs = 50, verbose = TRUE)
```

> PCA also takes a while on this many cells. Set `verbose = TRUE` so you can see
> it working (it prints the top genes for each PC as it runs).

**Your task:** Complete **Step 2a** in the template.

The `ElbowPlot` (written for you; it prints to the Plots pane — the `ggsave`
that would save `gata1_elbow.png` is commented out) shows how much variance each
PC explains; where the curve flattens is approximately where additional PCs
become noise. For this dataset the elbow is early, so we keep 10 PCs (20 also
works). Too few PCs merge distinct states; too many add noise. The aim is for
your conclusions to hold regardless of the exact choice.

```r
n_pcs <- 10
```

**Your task:** Complete **Step 2b** in the template.

---

## Step 5 — UMAP and graph-based clustering

- `RunUMAP(dims = 1:n_pcs)` projects the `n_pcs`-dimensional PCA space into 2D
  for visualization. UMAP is for visualization, not measurement: distances
  between well-separated groups are not quantitative.
- `FindNeighbors` followed by `FindClusters(resolution = 0.4)` performs the
  clustering on the PCA graph. Higher resolution produces more, smaller clusters.
  0.4 is a moderate choice; test a few values and select the one that yields
  biologically sensible groups.

```r
combined <- RunUMAP(combined, dims = 1:n_pcs)
combined <- FindNeighbors(combined, dims = 1:n_pcs)
combined <- FindClusters(combined, resolution = 0.4)
```

**Your task:** Complete **Steps 3a, 3b, and 3c** in the template.

> Cluster numbers are assigned arbitrarily: re-running `FindClusters` produces
> approximately the same groups, but the number assigned to each group may
> change.

---

## Step 6 — Visualize with `DimPlot`

`DimPlot` draws the UMAP colored by any metadata column. The pattern is always
the same:

```r
p <- DimPlot(combined, reduction = "umap", group.by = <a metadata column>)
p
save_dim(p, "some_name.png")
```

The bare `p` line prints the plot to the Plots pane. For the workshop, the
`save_dim()` helper is a **no-op** — its `ggsave` is commented out, so it just
views the plot. The filenames in the table below are what each call *would* save
if you re-enable the `ggsave` inside `save_dim`.

The first plot is provided as a worked example (colored by `seurat_clusters`).
You then repeat the pattern, changing only `group.by`:

| Template step | `group.by =` | Would save to |
|---|---|---|
| 3e | `"day"` | `gata1_umap_day.png` |
| 3f | `"construct"` | `gata1_umap_construct.png` |
| 3g | `"genotype"` | `gata1_umap_genotype.png` |
| 3i | `"Phase"` | `gata1_umap_ccPhase.png` |

The `Phase` plot (Step 3i) addresses a key question: are the clusters driven by
the cell cycle? If a cluster consists mostly of S/G2M cells, revisit the optional
regression in Step 3.

Two plots between these are written for you: the `S.Score` vs `G2M.Score` hex
scatter (Step 3h) and the `ApopScore1` UMAP (drawn with `FeaturePlot`, since
`ApopScore1` is a continuous score rather than a discrete group).

Finally, facet the day-colored UMAP into one panel per condition:

```r
p <- DimPlot(combined, reduction = "umap", group.by = "day",
             split.by = "construct_genotype", label = TRUE, pt.size = 0.3)
```

**Your task:** Complete **Steps 3e, 3f, 3g, 3i, and 3j** in the template. All the
`group.by` and `split.by` columns come from the metadata parsed in Lesson 01, so
each split is a single argument.

---

## Step 7 — Composition analysis

A UMAP can suggest that GATA1s shifts cells toward a lineage, but quantifying the
effect requires numbers. Cross-tabulate cluster × sample, then normalize the
table two ways, because each normalization answers a different question.

```r
ct <- table(Idents(combined), combined$sample)
```

**Your task:** Complete **Step 4a** in the template.

| Normalization | Code | Question answered |
|---|---|---|
| By cluster (rows sum to 1) | `prop.table(ct, margin = 1)` | Within this cluster, what fraction came from each sample? |
| By sample (columns sum to 1) | `prop.table(ct, margin = 2)` | Within this sample, what fraction of cells fell in each cluster? |

The by-sample version reveals proportion shifts across genotype, construct, and
day — for example, an erythroid cluster expanding in one condition.

**By cluster (margin = 1).** Build the table, convert it to a data frame, and
name the columns:

```r
frac_by_cluster <- prop.table(ct, margin = 1)
frac_df <- as.data.frame(frac_by_cluster)
colnames(frac_df) <- c("cluster", "sample", "fraction")
```

**Your task:** Complete **Steps 4b and 4c** in the template.

**By sample (margin = 2).** The same steps with a different margin and column
name:

```r
frac_by_sample <- prop.table(ct, margin = 2)
frac_df_sample <- as.data.frame(frac_by_sample)
colnames(frac_df_sample) <- c("cluster", "sample", "fraction_sample")
```

**Your task:** Complete the first three lines of **Step 4e** in the template.

The two `ggplot(... geom_tile() ...)` heatmaps (steelblue for by-cluster,
firebrick for by-sample) are written for you.

---

## Step 8 — Hierarchical clustering to order the heatmap

Plotting clusters and samples in arbitrary order obscures structure. We compute
Euclidean distances between rows and between columns, run `hclust`, and reorder
the heatmap axes by the resulting dendrogram so that similar samples are placed
together. This is the same principle as the dendrograms on a bulk RNA-seq
heatmap: the ordering is determined by the data.

The distance and `hclust` calls are written for you:

```r
m <- as.matrix(frac_by_cluster)
hc_col    <- hclust(dist(t(m), method = "euclidean"), method = "complete")
col_order <- hc_col$labels[hc_col$order]
```

Apply that ordering by converting the axes into factors whose level order
follows the clustering:

```r
frac_df$sample  <- factor(frac_df$sample,  levels = col_order)
frac_df$cluster <- factor(frac_df$cluster, levels = sort(unique(frac_df$cluster)))
```

**Your task:** Complete **Step 4d** (by-cluster) and the two `factor(...)` lines
of **Step 4f** (by-sample, using `row_order` for cluster and `col_order` for
sample). Both heatmaps print to the Plots pane (via `save_dim`, whose `ggsave` is
disabled); re-enable that `ggsave` to write `gata1_composition_by_cluster.png`
and `gata1_composition_by_sample.png`.

---

## Step 9 — Save for the annotation lesson (already written)

```r
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_clustered.rds"))
```

This clustered object is the input to the annotation lesson.

---

## Common pitfalls

- **Interpreting UMAP distances quantitatively.** Cluster membership is real; the
  2D geometry between clusters is not. Always support visual impressions with the
  composition tables.
- **Ignoring cell cycle.** If a cluster is entirely S/G2M cells, it may be a
  cell-cycle artifact rather than a cell type. Check the `Phase` UMAP and regress
  out `S.Score` and `G2M.Score` if needed.
- **Over-tuning resolution.** Increasing `FindClusters` resolution until you
  obtain the number of clusters you expected is circular reasoning. Choose a
  resolution, then validate with markers (the annotation lesson).
- **Reading the wrong margin.** `margin = 1` and `margin = 2` answer different
  questions. Confusing them flips the entire interpretation.
- **Skipping the regression.** If you do not regress out depth and mitochondrial
  fraction, the first PC may simply represent library size, and clusters will
  reflect technical variation.

## Check your understanding

1. Why log-normalize before PCA? What artifact does it remove?
2. What three columns does `CellCycleScoring` add, and what would you examine to
   decide whether cell cycle is a problem for your clusters?
3. You raise `FindClusters` resolution from 0.4 to 1.2. What happens to the
   cluster count, and how would you decide which resolution is correct?
4. A reviewer says an erythroid cluster expands with GATA1s. Which composition
   normalization (`margin = 1` or `margin = 2`) supports that claim, and why?
5. What does hierarchical clustering do to the heatmap axes, and why is that
   better than alphabetical order?
6. If the elbow plot showed variance still dropping steeply at PC 25, would you
   keep `n_pcs = 10`? What would you change?
