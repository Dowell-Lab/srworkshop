# GATA1 Lesson 03 — Multi-reference SingleR annotation and marker checks

**Template script (you fill this in):** `03_multiref_singler_annotation.R`
**Answer key:** `scripts_finished/03_multiref_singler_annotation.R`
**Input:** `OUT_DIR/gata1_combined_clustered.rds` (from Lesson 02)
**Output:** `OUT_DIR/gata1_combined_annotated_joined.rds` (the input to Lesson 04)

## How to use this worksheet

Open the template `03_multiref_singler_annotation.R` in RStudio next to this
worksheet. We work through the script one section at a time. For each step you
will see:

1. what the step does and why,
2. the code to write, and
3. a task cue indicating which `# YOUR CODE HERE:` block in the template to
   complete.

The boilerplate (loading libraries and data, the `run_singler` helper, the
`save_dim` helper, the marker gene list, and the side-by-side plot) is already
written for you. You complete the key line or lines in each step. Attempt each
one before consulting the answer key in `scripts_finished/`.

Each step provides the function name and its key arguments, so you assemble the
call rather than write it from scratch.

> Do not edit paths in the script. All paths live in
> `~/srworkshop/projectA/00_paths_and_setup.R`, and everything you create is
> written to `OUT_DIR`, not the read-only share.

## Overview

Lesson 02 produced clusters — but a cluster is just a numbered group. In this
lesson we give those groups biological *names* automatically with **SingleR**,
see why you might consult *several* reference atlases, and then verify the labels
against marker genes you trust. Automated labels are a hypothesis; markers are
the check. The overall flow is:

> JoinLayers → pull the data matrix → SingleR vs. a reference atlas →
> collapse rare labels → write metadata → UMAP by label → validate with markers →
> save

---

## Step 0 — Setup and load (already written)

There is nothing to complete here; read and run this block. The script sources
the shared paths file, loads the libraries (note the two new ones, `SingleR` and
`celldex`), and reads the clustered object saved by Lesson 02.

```r
source("~/srworkshop/projectA/00_paths_and_setup.R")
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(patchwork)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_clustered.rds"))
```

> If you just finished Lesson 02 in the same R session, `combined` is already in
> your environment. Otherwise reload it as above (or, if you never made it, use
> the pre-baked copy in `COOKING` — the commented `readRDS` line right below).

---

## Step 1 — One matrix for SingleR (`JoinLayers`)

After `merge()`, the RNA assay stores a separate `data` layer per sample
(`data.EuploidGATA1sD7`, `data.T21wtGATA1D9`, …). The script prints them with
`Layers(combined[["RNA"]])` so you can see the split. SingleR needs **one**
matrix of all cells, so we call `JoinLayers()` to stitch them into a unified
`data` layer, then pull it out with `LayerData(..., layer = "data")`. We feed the
log-normalized values, which is what SingleR's correlation scoring expects.

```r
combined_joined <- JoinLayers(combined)
sce_counts      <- LayerData(combined_joined, assay = "RNA", layer = "data")
```

**Your task:** Complete **Steps 1a and 1b** in the template.

---

## Step 2 — Why multiple references

`SingleR` labels each cell by correlating its profile against a reference of
known cell types ([Aran et al., *Nat Immunol* 2019](https://www.nature.com/articles/s41590-018-0276-y)).
The label you get depends on the reference you pick. `celldex` ships several,
and for iPSC-derived blood each has trade-offs:

| Reference | Strength for this dataset |
|---|---|
| `HumanPrimaryCellAtlasData()` | Broad human primary cells — good first pass |
| `BlueprintEncodeData()` | Broad stroma + immune |
| `MonacoImmuneData()` | Fine-grained immune subsets |
| `NovershternHematopoieticData()` | Classic hematopoietic hierarchy — blood-specific (our default here) |
| `DatabaseImmuneCellExpressionData()` | DICE immune populations |

`searchReferences()` lists more. The teaching point: **run a broad reference,
then a tissue-specific one, and compare.** Agreement between independent
references is far stronger evidence than one label alone. The template ships with
every option commented out and `NovershternHematopoieticData()` active; a
commented block in Step 4 lets you add a second reference and cross-tabulate the
two label sets.

> Heads up: loading a reference downloads it from the `celldex` cache the first
> time — it is **not** part of the read-only data share, so the first run can
> pause while it fetches. Let it finish.

The tiny `run_singler()` wrapper is written for you — it hands SingleR the test
matrix, the reference's `logcounts`, and the reference's main labels. You just
call it on your matrix:

```r
pred_broad <- run_singler(ref_broad, sce_counts)
table(pred_broad$labels)
```

**Your task:** Complete **Step 2a** in the template.

---

## Step 3 — Confidence, the long tail, and writing metadata

SingleR returns more than a label:

- **`pruned.labels`** — labels where low-confidence calls are set to `NA`. Prefer
  these for downstream analysis.
- **`delta.next`** — the score gap to the runner-up label. A small gap means the
  top two types were nearly tied, so trust that cell's call less.

With several references you also accumulate a long tail of labels assigned to
just a few cells, which clutters the UMAP. We keep labels with **> 150 cells**
and bucket the rest as `"other"` (`labels_other`) for a readable plot, while
retaining the full label in `SingleR_label`.

```r
pred_df <- as.data.frame(pred_broad) %>%
  add_count(labels, name = "label_n") %>%
  mutate(labels_other = case_when(label_n > 150 ~ labels, TRUE ~ "other"))
```

**Your task:** Complete **Step 3a** in the template.

We then write `SingleR_label`, `SingleR_labels_other`, `SingleR_pruned`, and
`SingleR_delta` back into the Seurat metadata (and copy `SingleR_labels_other`
onto `combined_joined`, since Lesson 04 reads the joined object). Three of those
assignments are written for you; you add the `SingleR_labels_other` column —
**the one CellChat groups by in Lesson 04.**

```r
combined$SingleR_labels_other <- pred_df$labels_other
```

**Your task:** Complete **Step 3b** in the template.

---

## Step 4 — Visualize the automated annotation

`DimPlot` draws the UMAP colored by any metadata column. The pattern is always
the same — only `group.by` changes:

```r
p <- DimPlot(combined, reduction = "umap", group.by = <a metadata column>)
p
save_dim(p, fn)
```

The first plot is provided as a worked example (colored by the full
`SingleR_label`). You then repeat the pattern for the pruned labels and the
readable `_other` labels:

| Template step | `group.by =` | Would save to |
|---|---|---|
| 5a (worked) | `"SingleR_label"` | `…gata1_umap_singler_label.png` |
| 5b | `"SingleR_pruned"` | `…gata1_umap_singler_pruned.png` |
| 5c | `"SingleR_labels_other"` | `…gata1_umap_singler_other.png` |

Comparing the three tells you a lot: where `SingleR_pruned` goes blank, SingleR
wasn't confident; where `_other` collapses a region into grey, that region was a
rag-bag of rare labels.

**Your task:** Complete **Steps 5b and 5c** in the template.

> Unlike the Lesson 01/02 `save_dim` (a no-op that only viewed the plot), here
> `save_dim` really writes a PNG to `OUT_DIR`. Every plot both prints to the
> Plots pane and lands as a file.

---

## Step 5 — Validate with canonical markers

Never trust an automated label you haven't sanity-checked. `FeaturePlot`s
overlay known genes on the UMAP:

- **`HBG2`, `HBG1`** — fetal hemoglobin → erythroid cells
- **`GATA1`** — the master erythroid/megakaryocyte TF, the gene of interest
- **`GFI1B`** — erythroid/megakaryocyte
- **`TFRC` (CD71)** — erythroid progenitors
- **`CD34`** — hematopoietic stem/progenitor cells
- **`PTPRC` (CD45)** — pan-leukocyte
- **`ITGA4`, `VAMP8`, `DIAPH3`** — additional markers used in the paper

The single-gene case (`agene <- "HBG2"`) is done for you as a worked example.
Then you paint the whole panel at once by handing `FeaturePlot` the full
`markers` vector (already defined for you):

```r
p_markers <- FeaturePlot(combined, features = markers, reduction = "umap",
                         pt.size = 0.3, ncol = 3)
```

**Your task:** Complete **Step 6a** in the template.

If the cluster SingleR called "erythroid" lights up for `HBG1`/`GATA1`/`TFRC`,
the label holds. If markers and labels disagree, believe the markers and
investigate. The script also saves a side-by-side of the labels next to `HBG1`
(written for you).

---

## Step 6 — Save for the CellChat lesson

The template ends with the two saves **commented out**, plus a `rm()` that frees
the large joined object. Lesson 04 reads `gata1_combined_annotated_joined.rds`,
so uncomment (at least) the second save before moving on:

```r
saveRDS(combined_joined, file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
```

This annotated object is the input to the cell-cell communication lesson.

---

## Common pitfalls

- **Skipping `JoinLayers`.** SingleR on a split assay errors or silently uses one
  layer. Join first.
- **Trusting one reference.** Different atlases name overlapping cell states
  differently. Compare at least two before committing.
- **Ignoring `pruned.labels`/`delta.next`.** The raw `labels` field always
  returns *something*, even for ambiguous cells. The confidence fields tell you
  which calls to distrust.
- **Believing labels over markers.** SingleR can only return cell types present
  in the reference. A novel or transitional state will be forced into the
  nearest reference label — markers catch this.
- **Forgetting to uncomment the save.** The saves ship commented out; if you skip
  them, Lesson 04 has nothing to read.

## Check your understanding

1. Why must you call `JoinLayers` before running SingleR on a merged object?
2. You annotate with a broad atlas and a hematopoietic atlas and they disagree
   for one cluster. How do you decide which to believe?
3. What do `pruned.labels` and `delta.next` each tell you, and how would you use
   them to filter cells before a proportion test?
4. The `> 150 cells` cutoff for `labels_other` is arbitrary. What is its purpose,
   and when might you change it?
5. A cluster is labeled "Monocytes" but is negative for `PTPRC` and positive for
   `GATA1`/`HBG1`. What is the most likely explanation, and what do you do?
