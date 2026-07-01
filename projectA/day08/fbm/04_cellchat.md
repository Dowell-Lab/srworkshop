# Lesson 4 — Cell-Cell Communication (CellChat)

**Script:** `scripts/04_cellchat.R`
**Dataset:** the **T21** (Trisomy 21) fetal bone marrow object (a matching **D21** control is available for comparison).
**Goal:** Infer which cell types are "talking" to each other through ligand-receptor signaling, and which pathways drive it.
**Vignette:** [CellChat tutorial](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html)

Clustering and annotation tell you *what* cells are present. CellChat asks the next question: *how do they coordinate?* It does this by checking, for every pair of cell types, whether a ligand made by one is matched by a receptor on the other.

The lesson is split because one step is genuinely slow:
- **Part 1** builds a CellChat object from a Seurat object (run once, ahead of time).
- **Part 2** analyzes a **pre-made** object — this is what you run in class.

---

## Part 1 — Building a CellChat object

### 2. `createCellChat(group.by = "cell.labels")`
CellChat needs to know each cell's type — that's the `group.by`. Communication is always computed *between cell-type groups*, so good labels (from Lesson 3) are a prerequisite.

### 3–4. The ligand-receptor database
`CellChatDB.human` is a curated catalog of signaling interactions, organized into categories (Secreted Signaling, ECM-Receptor, Cell-Cell Contact). `showDatabaseCategory()` and `glimpse()` let you inspect it. `subsetDB(..., "Secreted Signaling")` keeps just secreted-ligand interactions — a common memory-saving move. You attach the DB to the object and `subsetData()` restricts the expression data to genes in the DB.

### 5. Over-expressed genes & interactions
`identifyOverExpressedGenes` finds genes differentially expressed across your cell types; `identifyOverExpressedInteractions` keeps L-R pairs where at least one partner is over-expressed. This is what narrows thousands of possible interactions down to the ones with signal. `future::plan("multisession")` parallelizes it.

### 6. `computeCommunProb` (the slow step)
This models the actual communication probability for each L-R pair between each pair of cell types, using a law-of-mass-action-style model on mean expression. `type = "triMean"` uses a trimmed mean → fewer but stronger, higher-confidence interactions. **This is the step that's too slow for class**, so we precompute it and load the result in Part 2.

---

## Part 2 — Analyzing a pre-made object

We load the precomputed object with `readRDS()` (it's a `.rds` file).

### 2. `filterCommunication(min.cells = 10)`
Drops inferred communications that rely on cell groups with fewer than 10 cells — too few to trust statistically.

### 3. `computeCommunProbPathway`
Individual L-R pairs are aggregated up to the **signaling pathway** level (e.g. all CCL-family pairs → the "CCL" pathway). Pathways are the more interpretable unit of analysis.

### 4. `aggregateNet`
Counts/sums the interactions between every pair of cell types, producing the network we visualize.

### 5. Whole-dataset circle plot (`netVisual_circle`)
Nodes are cell types; edges are interactions; thicker edge = stronger communication. This is deliberately a mess — *everything* is plotted at once. The takeaway: you need a way to pick *which pathway* to focus on. That's the next step.

### 6. Centrality heatmaps (`netAnalysis_signalingRole_heatmap`)
After `netAnalysis_computeCentrality`, two heatmaps show, per pathway (rows) and cell type (columns), how strong the **outgoing** vs **incoming** signaling is (darker = stronger). **This is your pathway-picking tool** — start exploratory analysis from the darkest rows.

### 7. Per-pathway deep dive (example: CCL)
Why this matters biologically: **4 of the ~225 chromosome-21 genes are interferon-related**, and people with Trisomy 21 carry an extra chr-21 — so IFN/cytokine signaling (IFN-II, MHC-I/II, CXCL, CCL, CD40, TNF) is a natural focus. We examine **CCL** four ways:
- **a. `netVisual_aggregate(..., layout="circle")`** — the CCL-only communication network.
- **b. `netVisual_heatmap`** — sender × receiver strength for CCL.
- **c. `netAnalysis_contribution`** — which individual L-R pairs drive the CCL pathway.
- **d. `netVisual_bubble`** — significant L-R pairs from chosen senders (`sources.use`) to receivers (`targets.use`); `levels(cellChat@idents)` gives the index of each cell type so you can pick numbers.

### Comparing conditions
Repeat Part 2 with the **D21** object, then use CellChat's `mergeCellChat()` and comparison functions to contrast Trisomy 21 vs control signaling — the real scientific payoff.

---

## Common pitfalls
- Annotating poorly upstream → garbage communication (it's all driven by `cell.labels`).
- Reading the whole-dataset circle plot as a result rather than a "where do I look?" prompt.
- Forgetting that interaction *probability* is a model output, not a measured rate — treat it as a ranked hypothesis generator.

## Check your understanding
1. Why is `group.by` (cell labels) the single most important input to CellChat?
2. What do the outgoing/incoming centrality heatmaps let you decide?
3. Why is interferon-related signaling especially interesting in a Trisomy 21 dataset?
