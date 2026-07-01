# GATA1 Lesson 04 — Cell-cell communication with CellChat

**Script:** `scripts/gata1/04_cellchat.R`
**Input:** `OUT_DIR/gata1_combined_annotated.rds` (from Lesson 03)

Clustering tells you *which* cell types are present; CellChat asks *how they
talk to each other*. It infers signaling between cell groups from the expression
of known ligand-receptor pairs ([Jin et al., *Nat Commun* 2021](https://www.nature.com/articles/s41467-021-21246-9)).

> **Difference from the FBM track.** The fetal-bone-marrow lesson loads a
> *pre-made* CellChat `.rds` from the shared folder because the inference step is
> slow. There is no pre-made CellChat object for the GATA1 data, so **here you
> build one end to end from your own annotated object.** You get to see every
> step, and you save the result so you never recompute it.

## Step 1 — Choose the grouping label

CellChat groups cells by a cell-type label. We use the SingleR labels from
Lesson 03 — specifically the `SingleR_labels_other` version, which folds rare
labels into `"other"`. We then drop `other`/`NA` cells so CellChat isn't asked
to model groups too small to be meaningful. Good labels in → sensible signaling
out; garbage labels produce confident nonsense.

## Step 2 — Build the CellChat object

```r
cellChat <- createCellChat(object = combined, group.by = "cc_group", assay = "RNA")
cellChat@DB <- CellChatDB.human
```

`createCellChat` reads the log-normalized data and the group labels.
`CellChatDB.human` is the curated database of human ligand-receptor
interactions, organized into signaling pathways (secreted signaling,
cell-cell contact, ECM-receptor). `showDatabaseCategory()` shows the breakdown.

## Step 3 — Pre-processing

- `subsetData()` restricts to genes that appear in the database (everything else
  is irrelevant to signaling).
- `identifyOverExpressedGenes()` finds genes enriched in each cell group.
- `identifyOverExpressedInteractions()` keeps ligand-receptor pairs where both
  partners are over-expressed in the relevant groups.

We also set `future::plan("multisession")` so the heavy step runs in parallel.

## Step 4 — Infer the network (the slow step)

```r
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)
```

- `computeCommunProb` estimates, for every pair of cell groups and every
  ligand-receptor pair, the probability of communication. `type = "triMean"` is
  a robust summary of a group's expression. **This is the expensive call** — run
  it once, then `saveRDS` (the script does), and reload for plotting.
- `filterCommunication(min.cells = 10)` discards edges involving groups with too
  few cells.
- `computeCommunProbPathway` rolls individual L-R pairs up to pathway level.
- `aggregateNet` counts edges to build the summary network.

## Step 5 — Whole-dataset overview

Two circle plots side by side: nodes are cell types (sized by cell count), edges
show the **number** of interactions and the **strength** of interactions between
each pair. This is your bird's-eye view of who signals to whom.

## Step 6 — Which pathways carry the signal

`netAnalysis_computeCentrality` + `netAnalysis_signalingRole_heatmap` show, per
pathway, which cell types are dominant **senders** (outgoing) and **receivers**
(incoming). This is how you spot, for example, a progenitor population driving a
particular pathway.

## Step 7 — Focus on one pathway

The pathways available depend on the cell types your run recovered, so the
script **prints `cellChat@netP$pathways` and visualizes the first one** rather
than hard-coding a pathway that might be absent. Swap in a pathway of interest.
For your chosen pathway you get a circle plot, a heatmap, and a
`netAnalysis_contribution` bar chart showing which L-R pairs drive it.

## Step 8 — Comparing conditions (exercise)

The most interesting GATA1 question is *how signaling differs between
constructs or genotypes*. The recipe: subset the annotated object by `construct`
(or `genotype`), build a CellChat object for each subset (Steps 2–4), then use
`mergeCellChat()` and CellChat's comparison functions. Left as an exercise so
you practice the full build twice.

## Common pitfalls

- **Feeding raw counts.** CellChat expects log-normalized data. `JoinLayers`
  first, then it reads the `data` layer.
- **Too many tiny groups.** Rare labels create unstable probabilities. Fold them
  into `other` and drop them (we do), and use `min.cells`.
- **Recomputing `computeCommunProb` every time.** It's the slow step — save the
  object and reload.
- **Hard-coding a pathway name.** If that pathway isn't in your result, the plot
  errors. Always check `cellChat@netP$pathways` first.
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
