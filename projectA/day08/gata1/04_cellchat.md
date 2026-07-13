# GATA1 Lesson 04 — Cell-cell communication with CellChat

**Template script (you fill this in):** `04_cellchat.R`
**Answer key:** `scripts_finished/04_cellchat.R`
**Input:** `OUT_DIR/gata1_combined_annotated_joined.rds` (from Lesson 03)
**Output:** `OUT_DIR/gata1_cellchat.rds` (+ signaling plots in `OUT_DIR`)

## How to use this worksheet

Open the template `04_cellchat.R` in RStudio next to this worksheet. We work
through the script one section at a time. For each step you will see:

1. what the step does and why,
2. the code to write, and
3. a task cue indicating which `# YOUR CODE HERE:` block in the template to
   complete.

The boilerplate (loading libraries and data, building the `data.input` matrix and
`meta` table, saving the object, and all the plotting — the overview circles, the
signaling-role heatmap, and the per-pathway loop) is already written for you. You
complete the key line or lines in each step. Attempt each one before consulting
the answer key in `scripts_finished/`.

Each step provides the function name and its key arguments, so you assemble the
call rather than write it from scratch.

> Do not edit paths in the script. All paths live in `00_paths_and_setup.R`, and
> everything you create is written to `OUT_DIR`, not the read-only share.

## Overview

Clustering (Lesson 02) tells you *which* cell types are present; annotation
(Lesson 03) *named* them; CellChat now asks *how they talk to each other*. It
infers signaling between cell groups from the expression of known
ligand-receptor pairs ([Jin et al., *Nat Commun* 2021](https://www.nature.com/articles/s41467-021-21246-9)).
The overall flow is:

> pick a cell-type label → build the CellChat object → attach the L-R database →
> find over-expressed genes / interactions → computeCommunProb (slow) →
> aggregate the network → overview circles → signaling-role heatmap →
> per-pathway plots

> **Difference from the FBM track.** The fetal-bone-marrow lesson loads a
> *pre-made* CellChat `.rds` from the shared folder because the inference step is
> slow. There is no pre-made CellChat object for the GATA1 data, so **here you
> build one end to end from your own annotated object.** You get to see every
> step, and you save the result so you never recompute it.

---

## Step 0 — Setup and load (already written)

There is nothing to complete here; read and run this block. The script sources
the shared paths file, loads the libraries (note the new one, `CellChat`), reads
the annotated *joined* object from Lesson 03, and runs `gc()` to free memory.

```r
source("~/srworkshop/projectA/00_paths_and_setup.R")
library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)

combined <- readRDS(file.path(OUT_DIR, "gata1_combined_annotated_joined.rds"))
combined   # look at how many cells you have now
gc()       # garbage collector: trims the script's memory footprint
```

> If you never made the Lesson 03 object, use the pre-baked copy in `COOKING` —
> the commented `readRDS` line right below the first one.

---

## Step 1 — Choose the grouping label

CellChat groups cells by a cell-type label. We use the SingleR labels from
Lesson 03 — specifically the `SingleR_labels_other` version, which folds rare
labels into `"other"`. We copy it into a `cc_group` column, then drop
`other`/`NA` cells so CellChat isn't asked to model groups too small to be
meaningful. Good labels in → sensible signaling out; garbage labels produce
confident nonsense.

```r
combined$cc_group <- combined$SingleR_labels_other
combined <- subset(combined, subset = !is.na(cc_group) & cc_group != "other")
```

**Your task:** Complete **Steps 1a and 1b** in the template.

The `data.input` (the log-normalized matrix) and the `meta` data frame that pair
each cell barcode with its `cc_group` are built for you right after.

---

## Step 2 — Build the CellChat object

`createCellChat` reads the log-normalized data and the group labels.
`CellChatDB.human` is the curated database of human ligand-receptor
interactions, organized into signaling pathways (secreted signaling,
cell-cell contact, ECM-receptor). `showDatabaseCategory()` (written for you)
shows the breakdown; you attach the database to the object's `@DB` slot.

```r
cellChat <- createCellChat(object = data.input, meta = meta, group.by = "cc_group")

CellChatDB <- CellChatDB.human
cellChat@DB <- CellChatDB
```

**Your task:** Complete **Steps 2a and 2b** in the template.

---

## Step 3 — Pre-processing

Three calls trim the problem down to genes and pairs that could actually signal:

- `subsetData()` restricts to genes that appear in the database (everything else
  is irrelevant to signaling).
- `identifyOverExpressedGenes()` finds genes enriched in each cell group.
- `identifyOverExpressedInteractions()` keeps ligand-receptor pairs where both
  partners are over-expressed in the relevant groups.

```r
cellChat <- subsetData(cellChat)
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
```

**Your task:** Complete **Steps 3a, 3b, and 3c** in the template.

---

## Step 4 — Infer the network (the slow step)

```r
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)
```

- `computeCommunProb` estimates, for every pair of cell groups and every
  ligand-receptor pair, the probability of communication. `type = "triMean"` is
  a robust summary of a group's expression.
- `filterCommunication(min.cells = 10)` discards edges involving groups with too
  few cells.
- `computeCommunProbPathway` rolls individual L-R pairs up to pathway level.
- `aggregateNet` counts edges to build the summary network.

> Heads up: `computeCommunProb` (Step 4a) is **the expensive call** — on a small
> server it can take a while. Run it once. Right after `aggregateNet`, the
> script `saveRDS`-es the object (written for you) so you never recompute it —
> reload that `.rds` whenever you just want to re-plot.

**Your task:** Complete **Steps 4a, 4b, 4c, and 4d** in the template.

---

## Step 5 — Whole-dataset overview (already written)

Two circle plots side by side: nodes are cell types (sized by cell count), edges
show the **number** of interactions and the **strength** of interactions between
each pair. This is your bird's-eye view of who signals to whom. Because CellChat
draws these with base-R graphics, the block opens a `png()` device, draws, and
closes it with `dev.off()` — written for you; it writes
`gata1_cellchat_overview_circle.png` to `OUT_DIR`.

---

## Step 6 — Which pathways carry the signal (already written)

`netAnalysis_computeCentrality` + `netAnalysis_signalingRole_heatmap` show, per
pathway, which cell types are dominant **senders** (outgoing) and **receivers**
(incoming). This is how you spot, for example, a progenitor population driving a
particular pathway. The two heatmaps are drawn side by side into
`gata1_cellchat_signaling_roles.png` (written for you).

---

## Step 7 — Focus on one pathway

The pathways available depend on the cell types your run recovered, so the script
**reads `cellChat@netP$pathways`** rather than hard-coding a pathway that might be
absent:

```r
available_pathways <- cellChat@netP$pathways
print(available_pathways)
```

**Your task:** Complete **Step 7a** in the template.

The loop that follows (written for you) walks every recovered pathway and saves,
for each, a circle plot and a `netAnalysis_contribution` bar chart showing which
L-R pairs drive it. To zoom in on one pathway of interest, swap the loop for a
single `pathways.show <- available_pathways[i]` of your choice.

---

## Step 8 — Comparing conditions (exercise)

The most interesting GATA1 question is *how signaling differs between
constructs or genotypes*. The recipe: subset the annotated object by `construct`
(or `genotype`), build a CellChat object for each subset (Steps 1–4), then use
`mergeCellChat()` and CellChat's comparison functions. Left as an exercise (the
commented block at the bottom of the script) so you practice the full build
twice.

---

## Common pitfalls

- **Feeding raw counts.** CellChat expects log-normalized data. We read the
  joined object from Lesson 03 and pull its `data` layer.
- **Too many tiny groups.** Rare labels create unstable probabilities. Fold them
  into `other` and drop them (we do), and use `min.cells`.
- **Recomputing `computeCommunProb` every time.** It's the slow step — save the
  object and reload.
- **Hard-coding a pathway name.** If that pathway isn't in your result, the plot
  errors. Always read `cellChat@netP$pathways` first.
- **Over-reading absolute probabilities.** CellChat scores are comparative, not
  physical rates. Use them to rank and compare, not as literal fluxes.

## Check your understanding

1. Why does the GATA1 track build CellChat from scratch while the FBM track
   loads a pre-made object?
2. What does `computeCommunProb` actually estimate, and why is it the step you
   save to disk?
3. Why fold rare SingleR labels into `other` before running CellChat?
4. In the signaling-role heatmap, what distinguishes an "outgoing" from an
   "incoming" cell type for a pathway?
5. Outline the steps to test whether a pathway is stronger in GATA1s than in
   wtGATA1.
