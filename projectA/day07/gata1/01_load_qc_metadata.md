# GATA1 Lesson 01 — Loading data, QC, metadata-from-names, and module scores

**Script:** `scripts/gata1/day07/01_load_qc_metadata.R`
**Dataset:** GSE271399 — Takasaki et al., *Stem Cell Reports* 2025

This is the dataset you run from start to finish. It is a clean, well-designed
study, which makes it a good place to learn the habits that keep single-cell
analysis honest. The companion fetal-bone-marrow track (in `scripts/fbm/`) uses
a messier real-world dataset on purpose — we contrast the two in the
"metadata pitfalls" note.

## The experiment in one paragraph

The authors took isogenic induced pluripotent stem cell (iPSC) lines and varied
three things in a factorial design: **genotype** (Euploid vs. Trisomy 21),
**construct** (full-length *wtGATA1* vs. the truncated *GATA1s* isoform), and
**differentiation day** (D7, D9, D11). That is 2 × 2 × 3 = **12 conditions**.
They flow-purified hematopoietic progenitors and ran single-cell RNA-seq to ask
how trisomy 21 and GATA1s — individually and together — reshape blood
development ([GEO GSE271399](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271399),
[PubMed PMID 40680731](https://pubmed.ncbi.nlm.nih.gov/40680731/)).

GATA1s matters because the short isoform is the same lesion found in transient
abnormal myelopoiesis and myeloid leukemia of Down syndrome ([Wechsler et al.,
*Nat Genet* 2002](https://www.nature.com/articles/ng955)).

## Step 1 — Read each sample from the shared folder

The files already live on the cluster at the path stored in `GATA1_DIR`
(set in `scripts/00_paths_and_setup.R`). **You do not download anything.**

Each sample is three gzipped files in 10x layout:

- `..._matrix.mtx.gz` — the sparse counts matrix (genes × cells)
- `..._barcodes.tsv.gz` — one cell barcode per line
- `..._features.tsv.gz` — one gene per line; column 2 is the gene symbol

`ReadMtx()` glues those three pieces back into a matrix. We loop over
`GATA1_SAMPLES` so the 12 names live in exactly one place. Each matrix becomes a
`Seurat` object, and we stamp the sample name onto every cell with
`obj$sample <- smpl`. That single line is what makes the next steps possible.

`merge()` stacks all 12 objects into `combinedori`. This is a plain merge with
**no batch correction** — appropriate here because the lines are isogenic and
processed together. (The FBM track shows when you *do* need integration.)

## Step 2 — Mitochondrial percentage

`PercentageFeatureSet(obj, pattern = "^MT-")` computes, per cell, the fraction
of reads coming from mitochondrial genes (human mito symbols start with `MT-`).
A dying or membrane-compromised cell leaks cytoplasmic mRNA but keeps its
mitochondria, so a high `percent.mt` is a classic low-quality flag
([Ilicic et al., *Genome Biol* 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0888-1)).

## Step 3 — QC violin plots

We save one violin plot per sample showing three QC metrics:

- `nFeature_RNA` — number of distinct genes detected per cell
- `nCount_RNA` — total UMIs (transcripts) per cell
- `percent.mt` — mitochondrial fraction

Look at the spread *per sample* before choosing thresholds. The plots go to
`OUT_DIR/gata1_qc_violin_plots/`.

## Step 4 — Filter cells

```r
nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15
```

- `> 200` genes: drop empty droplets and debris.
- `< 6000` genes: drop probable doublets (two cells in one droplet read as one
  unusually rich "cell").
- `< 15%` mito: drop dying cells.

These are reasonable starting values, not laws of nature. The script writes
`gata1_qc_summary.csv` so you can see exactly how many cells each sample kept.

## Step 5 — Metadata FROM the sample names (the key habit)

The sample name `T21wtGATA1D9` already contains the entire design. Rather than
maintain a separate spreadsheet that someone has to keep in sync, we **parse**
the design straight out of the name with code:

```r
combined$genotype  <- ifelse(grepl("^Euploid", combined$sample), "Euploid", "T21")
combined$construct <- case_when(grepl("wtGATA1", ...) ~ "wtGATA1",
                                grepl("GATA1s",  ...) ~ "GATA1s")
combined$day       <- case_when(grepl("D7$", ...) ~ "D7", ...)
```

Because the name is the single source of truth, a sample **cannot** end up with
two contradictory labels. That is exactly the failure mode you will see in the
fetal-bone-marrow data, where one fetus is annotated with two impossible ages.
The `stopifnot()` checks at the end fail loudly if any cell didn't get a label —
catching typos in the sample list before they become silent bugs.

We also build combined factors (`day_construct`, `construct_genotype`, etc.) so
later plots can split panels by any slice of the design.

## Step 6 — Stress & apoptosis module scores

`AddModuleScore()` reduces a gene *set* to a single number per cell: the average
expression of the set minus a matched random-gene background
([Tirosh et al., *Science* 2016](https://www.science.org/doi/10.1126/science.aad0501)).
We score two sets:

- **stress** — immediate-early/heat-shock genes (`FOS`, `JUN`, `HSPA1A`, …)
  that spike during tissue dissociation.
- **apoptosis** — `BAX`, `CASP3`, `FAS`, … marking cells on the way out.

Comparing these across samples (and against `percent.mt`) confirms whether a
suspect sample such as `T21GATA1sD11` is genuinely stressed/dying rather than
biologically novel. This is a quantitative way to justify down-weighting or
excluding a sample.

The script saves `gata1_combined_qc.rds` for Lesson 02.

## Common pitfalls

- **Editing paths in the script.** Don't. All paths live in
  `00_paths_and_setup.R`; the scripts `source("../00_paths_and_setup.R")`.
- **Forgetting `obj$sample`.** Without it, every downstream `group.by`/`split.by`
  and the metadata parsing in Step 5 falls apart.
- **Treating QC thresholds as universal.** A 15% mito cap suits these cultured
  progenitors; primary tissue or other cell types may need different values.
  Always look at the violins first.
- **Reading a module score as absolute.** It is relative to a random background,
  so compare scores *between groups*, not against zero.
- **Writing to the shared folder.** It is read-only by design. Everything you
  create goes to `OUT_DIR`.

## Check your understanding

1. Why parse genotype/construct/day from the sample name instead of joining a
   metadata table? What error does this design make impossible?
2. A cell has `percent.mt = 40%` and only 150 genes detected. Keep or drop, and
   why?
3. What does a module score subtract, and why does that make it more trustworthy
   than a raw gene-set average?
4. You see `nFeature_RNA` up near 7000 for a cluster of cells. What artifact
   should you suspect, and which filter addresses it?
5. Why is a plain `merge()` (no integration) defensible for this isogenic
   dataset but risky for samples from different donors and batches?
