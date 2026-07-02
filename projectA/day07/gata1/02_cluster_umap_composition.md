# GATA1 Lesson 02 — Normalization, PCA, clustering, UMAP, and composition

**Script:** `scripts/gata1/02_cluster_umap_composition.R`
**Input:** `OUT_DIR/gata1_combined_qc.rds` (from Lesson 01)

Lesson 01 gave us a filtered object with clean metadata. Now we turn raw counts
into a map of cell states (clusters on a UMAP) and then *quantify* how the
experimental conditions are distributed across that map.

## Step 1 — Normalize → variable genes → scale

- **`NormalizeData` (LogNormalize, scale factor 1e4).** Cells differ in total
  counts purely for technical reasons. We divide each count by the cell's total,
  scale to 10,000, and log-transform so that "expression" is comparable across
  cells.
- **`FindVariableFeatures` (vst, 3000 genes).** Most genes are uninformative
  housekeepers. We keep the 3000 most variable genes, where the biology lives.
- **`ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt"))`.** Centers and
  scales each gene to mean 0, variance 1, and *regresses out* sequencing depth
  and mitochondrial fraction so those technical drivers don't masquerade as
  biological structure in the PCA.

## Step 2 — PCA and choosing the number of PCs

`RunPCA` compresses 3000 genes into a few dozen principal components. The
`ElbowPlot` (saved as `gata1_elbow.png`) shows how much variance each PC
explains; where the curve flattens is roughly where added PCs become noise. For
this dataset the elbow is early, so we keep **`n_pcs <- 10`**. Picking too few
PCs merges distinct states; too many adds noise. Re-examine the elbow if you
adapt this to other data.

## Step 3 — UMAP and graph-based clustering

- `RunUMAP(dims = 1:10)` projects the 10-PC space into 2D for visualization.
  **UMAP is for looking, not for measuring** — distances between far-apart blobs
  aren't quantitative.
- `FindNeighbors` + `FindClusters(resolution = 0.4)` do the actual clustering on
  the PCA graph. Higher resolution → more, smaller clusters. 0.4 is a moderate
  choice; try a couple of values and see which gives biologically sensible
  groups.

We save UMAPs colored by `seurat_clusters`, `day`, `construct`, `genotype`, and
a faceted `day` split by `construct_genotype`. Because the metadata was parsed
in Lesson 01, these splits are one argument each.

## Step 4 — Composition analysis (the quantitative payoff)

A UMAP can *suggest* that GATA1s shifts cells toward a lineage, but you need
numbers to claim it. We cross-tabulate cluster × sample with
`table(Idents(combined), combined$sample)` and then normalize two ways — each
answers a different question:

| Normalization | Code | Question answered |
|---|---|---|
| **By cluster** (rows sum to 1) | `prop.table(ct, margin = 1)` | "Within this cluster, what fraction came from each sample?" |
| **By sample** (columns sum to 1) | `prop.table(ct, margin = 2)` | "Within this sample, what fraction of cells fell in each cluster?" |

The **by-sample** version is the one that reveals proportion shifts across
genotype/construct/day — e.g. an erythroid cluster expanding in one condition.

### Hierarchical clustering to order the heatmap

Plotting clusters/samples in arbitrary order hides structure. We compute
Euclidean distances between rows and between columns, run `hclust`, and reorder
the heatmap axes by the dendrogram so similar samples sit together:

```r
hc_col    <- hclust(dist(t(mat), method = "euclidean"), method = "complete")
col_order <- hc_col$labels[hc_col$order]
frac_df$sample <- factor(frac_df$sample, levels = col_order)
```

This is the same idea as the dendrograms on a bulk RNA-seq heatmap: let the data
decide the ordering. Two heatmaps are saved —
`gata1_composition_by_cluster.png` (steelblue) and
`gata1_composition_by_sample.png` (firebrick).

The script saves `gata1_combined_clustered.rds` for Lesson 03.

## Common pitfalls

- **Interpreting UMAP distances quantitatively.** Cluster *membership* is real;
  the 2D geometry between clusters is not. Always back visual impressions with
  the composition tables.
- **Over-tuning resolution.** Cranking `FindClusters` resolution until you get
  the number of clusters you "expected" is circular. Choose a resolution, then
  validate with markers (Lesson 03).
- **Reading the wrong margin.** `margin = 1` vs. `margin = 2` answer different
  questions. Mixing them up flips the entire interpretation.
- **Skipping the regression.** If you don't regress out depth/mito, the first PC
  can just be "library size," and clusters reflect technical variation.

## Check your understanding

1. Why log-normalize before PCA? What artifact does it remove?
2. You raise `FindClusters` resolution from 0.4 to 1.2. What happens to cluster
   count, and how would you decide which resolution is "right"?
3. A reviewer says "an erythroid cluster expands with GATA1s." Which composition
   normalization (`margin = 1` or `margin = 2`) supports that claim, and why?
4. What is hierarchical clustering doing to the heatmap axes, and why is that
   better than alphabetical order?
5. If the elbow plot showed variance still dropping steeply at PC 25, would you
   keep `n_pcs = 10`? What would you change?
