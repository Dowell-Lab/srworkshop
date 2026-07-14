# Lesson 2 — Marker Genes & Doublet Detection

**Script:** `scripts/02_markers_doublets.R`
**Input:** `integrated_clustered.RData` from Lesson 1 (T21 + D21).
**Goal:** Find the genes that define each cluster, compare conditions, and remove doublets.

Once you have clusters, two questions follow immediately: *what are these clusters?* (markers) and *are any of these "cells" actually artifacts?* (doublets).

---

> **Memory tip — `gc()`.** R doesn't hand memory back to your computer the moment
> you delete or overwrite a large object; calling `gc()` (R's *garbage
> collector*) forces that cleanup and prints how much memory is now in use.
> Single-cell objects are big, so call `gc()` at the **heavy transitions** —
> right after `rm()`-ing an object you're done with, after a `merge()` /
> integration, or just before a memory-hungry step (SingleR, CellChat, monocle3).
> It never changes your results and is safe to run anytime; you just don't need it
> after every line.

---

## Part A — Differential expression / marker genes

### Why `JoinLayers` first?
During integration Seurat keeps each sample's counts in a separate **layer**. Marker testing needs to see all cells in one matrix, so `JoinLayers()` collapses them back together. Skip this and `FindAllMarkers` may error or test the wrong thing.

### 1. Markers across all clusters (`FindAllMarkers`)
For each cluster, this tests every gene: *is it higher in this cluster than in all other cells combined?* The output table has one row per significant gene per cluster, with key columns:
- **`avg_log2FC`** — effect size (how much higher, log2 scale).
- **`pct.1` / `pct.2`** — fraction of cells expressing the gene inside / outside the cluster.
- **`p_val_adj`** — multiple-testing-corrected p-value; trust this, not raw `p_val`.

We rank two ways:
- **Top by `avg_log2FC`** → the genes that most strongly *define* the cluster.
- **Top by `pct.1`** → clean, broadly-detected markers (good for FeaturePlots / labeling).

### 2. Markers between two specific clusters (`FindMarkers`)
`ident.1` vs `ident.2` is a direct head-to-head. Positive `avg_log2FC` = up in `ident.1`. Use this when two clusters look similar and you want to know what separates them.

### 3. Differential expression by *condition* (metadata)
This is the part most relevant to the T21 vs D21 question. `SetIdent(value = "T21.status")` re-labels every cell by condition instead of cluster, so `FindAllMarkers` now finds **genes that differ between Trisomy 21 and Disomic cells**.

> **Crucial nuance:** comparing all T21 vs all D21 confounds *cell-type composition* with *expression*. If T21 has more erythroid cells, erythroid genes look "up in T21" even if no gene changed within any cell type. The fix shown: **subset to one cluster first, then compare T21 vs D21 inside it.** That controls for composition.

### 4. Visualizing a marker (`FeaturePlot`)
`FeaturePlot` colors the UMAP by a single gene's expression. `SLC4A1` lights up erythroid cells; `FCGR3B` lights up neutrophils. Always `DefaultAssay(obj) <- "RNA"` first so you plot measured expression, not integrated/scaled values.

---

## Part B — Doublet detection with DoubletFinder

### What is a doublet?
Two cells captured in one droplet share a barcode, producing a single "cell" whose profile is a blend of two real types. Doublets create fake intermediate populations and spurious co-expression — they must be flagged and usually removed.

### How DoubletFinder works
It **spikes in artificial doublets** (by averaging random pairs of real cells), projects everything into PCA space, and asks: *which real cells have the most artificial doublets as neighbors?* Those are the suspects.

### The parameters
- **`pN`** — proportion of artificial doublets to generate (0.25 is the standard default; results are insensitive to it).
- **`pK`** — neighborhood size. This one **matters**; ideally tune it with `paramSweep` + `find.pK` rather than the fixed `0.09` used here for speed.
- **`nExp`** — how many doublets you *expect*. Driven by loading density; ~10% is a reasonable rule of thumb but check your 10x loading concentration.
- **`PCs`** — which principal components to use.
- **`sct`** — set `TRUE` only if you normalized with SCTransform; we used log-normalization, so `FALSE`.

### Reading the result
The classification lands in a metadata column like `DF.classifications_0.25_0.09_NNN`. We grab its name with `grepl` (the suffix changes per run), then color the UMAP by Singlet/Doublet. Doublets often sit *between* clusters — exactly where blended profiles would land. The final commented line shows how to subset down to singlets.

---

## Common pitfalls
- Forgetting `JoinLayers` before marker testing.
- Interpreting whole-dataset T21-vs-D21 DE without controlling for composition.
- Using the default `pK` blindly — tune it for real analyses.
- Setting `sct = TRUE` when you didn't actually run SCTransform.

## Check your understanding
1. A gene appears "up in T21" across the whole dataset but is flat within every cluster. What happened?
2. Why do doublets tend to appear between clusters on a UMAP?
3. When would you prefer ranking markers by `pct.1` over `avg_log2FC`?

## Output
`all_clusters_diff_genes.RData` (marker table) saved to `OUT_DIR`.
