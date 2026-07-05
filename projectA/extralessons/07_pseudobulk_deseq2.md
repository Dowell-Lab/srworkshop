# Advanced 07 — Pseudobulk differential expression with DESeq2

**Script:** `scripts/advanced/07_pseudobulk_deseq2.R`
**Dataset:** GATA1 (GSE271399), annotated object from `gata1/03`

---

## Why this lesson exists

When you want to compare **conditions** in single-cell data — say GATA1s vs wtGATA1 — the tempting move is to run `FindMarkers` across the two conditions and read off the p-values. **This is statistically wrong**, and it is one of the most common mistakes in the field.

The problem: thousands of cells from **one sample** are *not* thousands of independent biological replicates. They share that sample's donor, its handling, its batch. Treating each cell as an independent observation massively inflates significance — you get thousands of "significant" genes that are really just one sample looking different from another.

**Pseudobulk** fixes this by restoring the true unit of replication: the sample.

> Squair et al., *Nature Communications* 2021 showed single-cell DE methods that ignore sample structure have badly inflated false-positive rates, and that pseudobulk methods (DESeq2, edgeR) are far better calibrated. https://www.nature.com/articles/s41467-021-25960-2

The GATA1 dataset is ideal for this because it has a **clean factorial design** (genotype × construct × day), giving you real biological replicates to aggregate.

---

## The core idea: collapse cells into pseudobulk

For each **(cell type × sample)** combination, we **sum the raw counts** across all its cells into a single profile — one "pseudobulk" sample. This turns a single-cell matrix into a small bulk-RNA-seq-like matrix with a handful of columns (your real replicates), which is exactly what DESeq2 was designed for.

`AggregateExpression(return.seurat = TRUE, group.by = ...)` does the summing. **Summed raw counts** — not averages, not normalized values — are what a negative-binomial model like DESeq2 expects.

---

## What the script does, step by step

1. **Load** the annotated GATA1 object and drop the `"other"` cell-type bucket so we test real cell types only.
2. **Aggregate** to pseudobulk: `AggregateExpression` summed counts per (cell type × sample × construct), then `round()` to guarantee integers.
3. **Filter genes** with edgeR's `filterByExpr(min.count = 10)` — removes genes too lowly expressed to be testable, which stabilizes the dispersion estimates.
4. **Build the DESeq2 dataset** with `DESeqDataSetFromMatrix(design = ~ construct)`, after `relevel`-ing so `wtGATA1` is the reference level (the denominator of the fold change).
5. **Run** `DESeq()`, wrapped in a `tryCatch` that falls back to `sfType = "poscounts"` — the size-factor estimator that survives genes with zeros, common in sparse pseudobulk.
6. **QC the samples** with `rlog(blind = TRUE)`: a PCA of pseudobulk samples plus a sample-to-sample distance heatmap (`pheatmap`, ward.D). Replicates should cluster by condition; an outlier here is a warning sign.
7. **Extract results**, plot a **p-value histogram** (a sanity check on calibration), and optionally refine with `fdrtool` and shrink fold changes with `lfcShrink(type = "ashr")`.
8. **MA plot** and a **merged results CSV** written to `OUT_DIR/advanced_pseudobulk_.../`.

---

## Reading the output

- **rlog PCA of pseudobulk samples** — the sanity check that matters most. Samples should separate by construct along a PC. If they separate by day or batch instead, add that covariate to the design (`~ day + construct`).
- **Distance heatmap** — replicates of the same condition should be each other's nearest neighbors.
- **p-value histogram** — should be roughly uniform with a spike near zero. A hump in the middle or a spike near 1 signals a model or filtering problem, not real signal.
- **MA plot** — after `ashr` shrinkage, low-count genes stop shouting large, noisy fold changes.

---

## Why the extra steps (edgeR, fdrtool, ashr)

- **edgeR `filterByExpr`** — a principled, well-tested gene filter; better than an arbitrary "expressed in N cells" rule.
- **`rlog`** — variance-stabilizing transform used *only for visualization/QC*, never for the test itself (the test uses raw counts).
- **`fdrtool`** — re-estimates the FDR from the empirical null when DESeq2's p-value distribution is mis-calibrated.
- **`lfcShrink(type = "ashr")`** — shrinks fold changes toward zero for low-information genes so your ranked gene list isn't dominated by noisy, low-count outliers. This matters when the shrunken log-fold-changes feed into GSEA (see lesson 08).

---

## Common pitfalls
- **Using `FindMarkers` across conditions and trusting the p-values.** Pseudo-replication inflates significance; this is the whole reason pseudobulk exists.
- **Aggregating averages or normalized data instead of summed raw counts.** DESeq2's model needs integer counts.
- **Too few replicates.** Pseudobulk trades cells for samples — if you only have one sample per condition, you have no replication and DESeq2 cannot estimate dispersion honestly.
- **Ignoring known covariates.** If day or donor drives the PCA, put it in the design (`~ day + construct`) or the construct effect will be confounded.
- **Reporting `rlog` values as the result.** `rlog` is for QC plots; the differential test is on counts.
- **Skipping the p-value histogram.** It is the fastest way to catch a broken model.

## Check your understanding
1. Why is running `FindMarkers` across two conditions a statistically flawed way to compare them?
2. What exactly is summed, and over what grouping, to create a pseudobulk sample?
3. Why must the aggregated values be summed raw counts rather than averaged or normalized expression?
4. What would you conclude if the rlog PCA separated pseudobulk samples by **day** rather than by **construct**, and how would you fix it?
5. What does `lfcShrink(type = "ashr")` accomplish, and why does it matter before feeding a gene ranking into GSEA?
