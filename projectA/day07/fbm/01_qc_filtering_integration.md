# Lesson 1 — QC, Filtering & Integration

**Script:** `scripts/01_qc_filtering_integration.R`
**Datasets:** `T21BM_male04` (Trisomy 21) and `D21_male35` (Disomic control), raw 10x CellRanger matrices on the shared AWS share.
**Goal:** Turn raw droplet counts into a clean, integrated, clustered Seurat object you can trust.

This is the foundation. Every later lesson assumes you understand what comes out of this one.

---

## The big picture

A single-cell experiment gives you a giant **genes × cells** count matrix. Before any biology, you have to answer three questions:

1. **Which "cells" are real?** (some droplets are empty or contain debris)
2. **Are two samples comparable?** (or is the difference between them just batch effect?)
3. **What groups of cells exist?** (clustering)

This lesson walks the chain: **read → QC → filter → normalize → integrate → cluster**.

---

## Step-by-step

### 1. Read the 10x matrices (`Read10X` → `CreateSeuratObject`)
`Read10X()` points at a CellRanger `filtered_feature_bc_matrix/` folder, which holds three files: barcodes (cell IDs), features (genes), and the sparse count matrix. `CreateSeuratObject()` wraps that matrix in the Seurat container that all downstream functions expect. The `project` argument sets `orig.ident`, the per-cell sample label we use later to color plots and split comparisons.

### 2. Ambient RNA correction with SoupX (shown, optional)
Some of the counts in every droplet come from free-floating mRNA in the cell suspension — the "soup" — not from the cell itself. **SoupX** estimates that contamination fraction and subtracts it. It needs the *raw* (unfiltered) matrix, the *filtered* matrix, and clustering, all of which `load10X()` reads from the CellRanger `outs/` folder. We show the call commented out so you know where it belongs; the rest of the lesson runs on the filtered matrices directly.

### 3. Add metadata
We tag every cell with `T21.status` ("T21" vs "D21") and `gender`. Metadata is how you later ask questions like "which genes differ between conditions?" — if it isn't in `@meta.data`, you can't group by it.

### 4. Percent mitochondrial reads (`PercentageFeatureSet`)
A cell that is dying or lysed leaks cytoplasmic mRNA but keeps its mitochondrial transcripts, so a **high `percent.mt` is a red flag**. Human mito genes start with `MT-` (the `^MT-` regex anchors to the start). The `VlnPlot` shows the three QC metrics that matter:
- **`nFeature_RNA`** — number of distinct genes per cell
- **`nCount_RNA`** — total transcripts per cell
- **`percent.mt`** — mitochondrial fraction

### 5. Filter (`subset`)
We drop cells outside sensible bounds:
- `nFeature_RNA > 200` removes empty droplets / debris (too few genes detected).
- `nFeature_RNA < 2500` removes likely **doublets** (two cells in one droplet → suspiciously many genes).
- `percent.mt < 5` removes dying cells.

> **Important teaching point:** these numbers are *not* universal. They were chosen by *looking at the violin plots*. Always set thresholds from your own QC distributions; the upper `nFeature` cut especially varies a lot by tissue.

We re-plot after filtering to confirm the extreme tails are gone.

### 6. Normalize → variable features → scale → PCA
This four-step chain is standard Seurat preprocessing:
- **`NormalizeData`** — log-normalizes so a cell with more total counts isn't automatically "higher" everywhere.
- **`FindVariableFeatures`** — picks the ~2000 most variable genes; these carry the biological signal, and using them speeds everything up.
- **`ScaleData`** — centers and scales each gene so PCA isn't dominated by a few highly expressed genes.
- **`RunPCA`** — linear dimensionality reduction; the PCs become the input to integration and clustering.

### 7. Integration with RPCA (`IntegrateLayers`)
If you just merged two samples and clustered, cells might separate **by sample** rather than **by cell type** — a batch effect. **RPCA (Reciprocal PCA) integration** aligns the shared populations across T21 and D21 so the same cell type from each sample overlaps. The result is a new reduction, `integrated.rpca`, that we cluster on instead of raw PCA.

> Why RPCA and not CCA? RPCA is faster and more conservative — it assumes the samples are fairly similar (true here: same tissue, differing only by trisomy), which avoids over-correcting away real biology.

### 8. Cluster + UMAP (`FindNeighbors` → `RunUMAP` → `FindClusters`)
- **`FindNeighbors`** builds a k-nearest-neighbor graph in the integrated space (`dims = 1:30` uses the first 30 integrated components).
- **`RunUMAP`** projects to 2D *for visualization only*. `n.neighbors` and `min.dist` control how tight/spread the layout looks — they don't change the underlying data.
- **`FindClusters`** runs Louvain community detection on the graph. **`resolution`** is the knob: higher → more, smaller clusters. `0.2` is intentionally low for a clean teaching view.

The `DimPlot`s let you eyeball success: clusters should reflect cell types, and T21/D21 cells should be *intermixed* within clusters (good integration), not split apart.

---

## Common pitfalls
- Clustering on `pca` instead of `integrated.rpca` → batch-driven clusters.
- Forgetting to re-plot QC after filtering → you never confirm the filter worked.
- Treating `resolution` as "correct/incorrect" — it's exploratory; try several.

## Check your understanding
1. Why might `nFeature_RNA` being *too high* indicate a doublet rather than a great cell?
2. If T21 and D21 cells formed two separate islands on the UMAP, what step would you suspect first?
3. What changes if you raise `resolution` from 0.2 to 1.0?

## Output
`integrated_clustered.RData` (in your writable `OUT_DIR`) → used by Lesson 2.
