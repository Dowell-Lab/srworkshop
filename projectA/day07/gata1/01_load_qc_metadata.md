# GATA1 Lesson 01 — Loading data, QC, metadata-from-names, and module scores

**Template script (you fill this in):** `01_load_qc_metadata.R`
**Answer key:** `scripts_finished/01_load_qc_metadata.R`
**Dataset:** GSE271399 — Takasaki et al., *Stem Cell Reports* 2025
**Output:** `OUT_DIR/gata1_combined_qc.rds` (the input to Lesson 02)

## How to use this worksheet

Open the template `01_load_qc_metadata.R` in RStudio next to this worksheet. We
work through the script one section at a time. For each step you will see:

1. what the step does and why,
2. the code to write, and
3. a task cue indicating which `# YOUR CODE HERE:` block in the template to
   complete.

The boilerplate (loading libraries, the QC violin-plot loop, the QC-summary CSV,
the module-score plots, saving the object) is already written for you. You
complete the key line or lines in each step. Attempt each one before consulting
the answer key in `scripts_finished/`.

Each step provides the function name and its key arguments, so you assemble the
call rather than write it from scratch.

> Do not edit paths in the script. All paths live in
> `~/srworkshop/projectA/00_paths_and_setup.R`, and everything you create is
> written to `OUT_DIR`, not the read-only share.

## Overview

This is the dataset you run from start to finish. It is a clean, well-designed
study, which makes it a good place to learn the habits that keep single-cell
analysis honest. The companion fetal-bone-marrow track uses a messier real-world
dataset on purpose — we contrast the two in the "metadata pitfalls" note. The
overall flow is:

> read 10x matrices → merge → %mito → QC filter → metadata-from-names →
> module scores → save

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

> The active `GATA1_SAMPLES` list in `00_paths_and_setup.R` may be a subset of
> the full 12 conditions (the day panel is trimmed to keep the exercise fast).
> Everything below works unchanged whichever samples are enabled, because the
> sample names live in exactly one place.

---

## Step 0 — Setup and libraries (already written)

There is nothing to complete here; read and run this block. The script sources
the shared paths file and loads the libraries.

```r
source("~/srworkshop/projectA/00_paths_and_setup.R")
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
```

> **On the AWS/Fiji cluster?** The `source(...)` line above is all you need — it
> already points `GATA1_DIR`, `OUT_DIR`, and `GATA1_SAMPLES` at the shared
> read-only data and a writable output folder. Skip the block below.
>
> **Working locally on your own machine?** Still run the `source(...)` line above
> (it defines the helpers), then **override the three directory variables** so
> they point at your local data and output folders. Run this right after
> `source(...)`:

```r
# Assumes your R working directory is the
# folder that contains data/ (set it with an RStudio Project or setwd()).
GATA1_DIR <- paste(getwd(), "data", sep = "/")     # your local 10x matrices

GATA1_SAMPLES <- c(
  "EuploidGATA1sD7",  "EuploidGATA1sD9",  "EuploidGATA1sD11",
  "EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
  "T21GATA1sD7",      "T21GATA1sD9",      "T21GATA1sD11",
  "T21wtGATA1D7",     "T21wtGATA1D9",     "T21wtGATA1D11"
)

OUT_DIR <- paste(getwd(), "outdir", sep = "/")     # where your results go
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)  # make it if missing
```

> **Memory tip — `gc()`.** R doesn't hand memory back to your computer the moment
> you delete or overwrite a large object; calling `gc()` (R's *garbage
> collector*) forces that cleanup and prints how much memory is now in use.
> Single-cell objects are big, so call `gc()` at the **heavy transitions** —
> right after `rm()`-ing an object you're done with, after a `merge()` /
> integration, or just before a memory-hungry step (SingleR, CellChat, monocle3).
> It never changes your results and is safe to run anytime; you just don't need it
> after every line.

---

## Step 1 — Read each sample and merge

The files already live on the cluster at the path stored in `GATA1_DIR`. **You
do not download anything.** Each sample is three gzipped files in 10x layout:

- `..._matrix.mtx.gz` — the sparse counts matrix (genes × cells)
- `..._barcodes.tsv.gz` — one cell barcode per line
- `..._features.tsv.gz` — one gene per line; column 2 is the gene symbol

`ReadMtx()` glues those three pieces back into a matrix. The template shows this
once as a **worked example** on a single sample (**Step 1a** — already written),
then loops over `GATA1_SAMPLES` so the sample names live in exactly one place.

> Heads up: that loop reads every sample and builds a Seurat object for each
> one, so it can take a little while to finish. Start it and be patient — it is
> done when your prompt returns.

**Build one object per sample.** Inside the loop, wrap each matrix in a Seurat
object and — critically — stamp the sample name onto every cell:

```r
obj <- CreateSeuratObject(counts = mat, project = smpl)
obj$sample <- smpl      # stamp the sample name onto every cell
```

That `obj$sample <- smpl` line is what makes every downstream `group.by` /
`split.by` and the metadata parsing in Step 5 possible.

**Your task:** Complete **Steps 1b and 1c** in the template.

**Merge.** `merge()` stacks all the objects into `combinedori`. This is a plain
merge with **no batch correction** — appropriate here because the lines are
isogenic and processed together. (The FBM track shows when you *do* need
integration.) We keep this pre-filter object so we can compare cell counts before
vs. after QC.

```r
combinedori <- merge(x = seurat_list[[1]], y = seurat_list[-1],
                     add.cell.ids = GATA1_SAMPLES)
```

**Your task:** Complete **Step 1d** in the template.

---

## Step 2 — Mitochondrial percentage

`PercentageFeatureSet(obj, pattern = "^MT-")` computes, per cell, the fraction of
reads coming from mitochondrial genes (human mito symbols start with `MT-`). A
dying or membrane-compromised cell leaks cytoplasmic mRNA but keeps its
mitochondria, so a high `percent.mt` is a classic low-quality flag
([Ilicic et al., *Genome Biol* 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0888-1)).

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)
```

**Your task:** Complete **Step 2a** in the template (the `lapply` wrapper is
provided; you write the one line inside it).

---

## Step 3 — QC violin plots (already written)

The script saves one violin plot per sample showing three QC metrics:

- `nFeature_RNA` — number of distinct genes detected per cell
- `nCount_RNA` — total UMIs (transcripts) per cell
- `percent.mt` — mitochondrial fraction

Look at the spread *per sample* before choosing thresholds. For the workshop the
violins simply print to the **Plots pane** (the `ggsave` lines are commented out;
uncomment them if you want PNG files in `OUT_DIR/gata1_qc_violin_plots/`). This
whole block is written for you — but **read the plots**, because Step 4 depends
on what you see. Two samples look off; part of the exercise is spotting which
and why.

---

## Step 4 — Filter cells

```r
subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
```

- `> 200` genes: drop empty droplets and debris.
- `< 6000` genes: drop probable doublets (two cells in one droplet read as one
  unusually rich "cell").
- `< 15%` mito: drop dying cells.

These are reasonable starting values, **not laws of nature** — you chose them by
looking at the violins in Step 3. The template wraps this in an `lapply`; you
write the `subset(...)` call inside. Everything after it (re-merging the filtered
objects and printing a QC summary table that reports how many cells each sample
kept) is written for you. The summary is printed to the console; the
`write.csv` line that would save it to `gata1_qc_summary.csv` is commented out —
uncomment it if you want the file.

**Your task:** Complete **Step 4a** in the template.

---

## Step 5 — Metadata FROM the sample names (the key habit)

The sample name `T21wtGATA1D11` already contains the entire design. Rather than
maintain a separate spreadsheet that someone has to keep in sync, we **parse**
the design straight out of the name with code. Because the name is the single
source of truth, a sample **cannot** end up with two contradictory labels — which
is exactly the failure mode you will see in the fetal-bone-marrow data, where one
fetus is annotated with two impossible ages.

```r
combined$genotype  <- ifelse(grepl("^Euploid", combined$sample), "Euploid", "T21")

combined$construct <- dplyr::case_when(
  grepl("wtGATA1", combined$sample) ~ "wtGATA1",
  grepl("GATA1s",  combined$sample) ~ "GATA1s",
  TRUE                              ~ NA_character_
)

combined$day <- dplyr::case_when(
  grepl("D7$",  combined$sample) ~ "D7",
  grepl("D9$",  combined$sample) ~ "D9",
  grepl("D11$", combined$sample) ~ "D11",
  TRUE                           ~ NA_character_
)
```

**Your task:** Complete **Steps 5a, 5b, and 5c** in the template.

The combined factors (`day_construct`, `construct_genotype`, …) that let later
plots split panels by any slice of the design, and the `stopifnot()` checks that
fail loudly if any cell didn't get a label, are written for you. Those checks
catch typos in the sample list before they become silent bugs.

---

## Step 6 — Stress & apoptosis module scores

`AddModuleScore()` reduces a gene *set* to a single number per cell: the average
expression of the set minus a matched random-gene background
([Tirosh et al., *Science* 2016](https://www.science.org/doi/10.1126/science.aad0501)).
We score two sets:

- **stress** — immediate-early/heat-shock genes (`FOS`, `JUN`, `HSPA1A`, …) that
  spike during tissue dissociation.
- **apoptosis** — `BAX`, `CASP3`, `FAS`, … marking cells on the way out.

`NormalizeData` is run first (written for you) so the scores read sensible values.
Then score each set:

```r
combined <- AddModuleScore(combined, features = list(stress_genes),    name = "StressScore")
combined <- AddModuleScore(combined, features = list(apoptosis_genes), name = "ApopScore")
```

Note the trailing `1` Seurat appends: the columns land as `StressScore1` and
`ApopScore1`. Comparing these across samples (and against `percent.mt`) confirms
whether a suspect sample is genuinely stressed/dying rather than biologically
novel — a quantitative way to justify down-weighting or excluding it.

**Your task:** Complete **Steps 6a and 6b** in the template. The `VlnPlot` that
compares the suspect sample to everyone else is written for you.

---

## Step 7 — Save for Lesson 02 (already written)

```r
saveRDS(combined, file.path(OUT_DIR, "gata1_combined_qc.rds"))
```

This filtered, annotated object is the input to the clustering lesson.

---

## Common pitfalls

- **Editing paths in the script.** Don't. All paths live in
  `00_paths_and_setup.R`; the scripts `source(...)` it.
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
