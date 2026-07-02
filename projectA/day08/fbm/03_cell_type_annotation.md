# Lesson 3 — Cell Type Annotation

**Script:** `scripts/03_cell_type_annotation.R`
**Dataset:** the pre-labeled **T21** (Trisomy 21) fetal bone marrow object.
**Goal:** Put biological names on clusters automatically, using three complementary strategies, and *measure* how well each works.

Clusters are just numbers until you name them. There are two broad philosophies for naming, and we cover both — plus a built-in-atlas shortcut.

| Approach | Needs a labeled reference? | Tool |
|---|---|---|
| A. Reference mapping | Yes (your own) | Seurat `TransferData` |
| B. Atlas scoring | Yes (built-in) | SingleR + celldex |
| C. Marker scoring | No | SC-Type |

### Why this dataset is great for teaching
The object already carries **true labels** at two resolutions — `broad_extfig7A_cell.labels` (coarse: B_lineage, erythroid, monocyte…) and `cell.labels` (fine sub-states). We **split it in half**, treat one half as a labeled "reference" and the other as an "unknown" query, then check predictions against the labels we hid. The `set.seed(42)` makes everyone's split identical.

> The script auto-detects the UMAP reduction name (`Xumap_` vs `umap`) and the label column, so it runs even if the object's internals differ slightly.

---

## Approach A — Seurat reference mapping

The logic is "label by neighbor": cells that look alike across datasets probably *are* alike.

1. **`FindTransferAnchors`** finds **anchors** — pairs of cells (one reference, one query) that are mutual nearest neighbors in a shared space. Anchors are the bridges across which labels travel.
2. **`TransferData`** uses those anchors to predict a label and a confidence score for every query cell. `refdata` is the column of reference labels you want to copy over.
3. **`AddMetaData`** (here a direct `$` assignment) attaches predictions to the query.
4. **Scoring:** because we know the truth, `predicted == true` gives an accuracy `table()`. We expect *high but not perfect* accuracy — it's half of the same dataset, but mapping still makes mistakes at cell-type boundaries.

**Teaching beat:** run it again with the *fine* `cell.labels` and watch accuracy drop. The closer two cell states are, the harder they are to tell apart — annotation quality depends on how distinct your clusters actually are.

---

## Approach B — SingleR + celldex

When you don't have your own labeled reference, **celldex** ships curated public atlases (e.g. `HumanPrimaryCellAtlasData()`, Monaco immune, Blueprint/ENCODE). **SingleR** correlates each query cell's expression against the reference's per-label expression profiles and assigns the best-matching label.

Key points:
- Feed SingleR the **log-normalized** data (`layer = "data"`), not scaled.
- `labels = hpca$label.main` uses coarse labels; `label.fine` gives finer ones.
- SingleR labels **per cell** (not per cluster), so it's robust to your clustering choices — a nice contrast to the cluster-level methods.
- Choose a reference that *matches your tissue*. A generic atlas will call bone marrow cells approximately, but a hematopoietic reference does better.

---

## Approach C — SC-Type (marker-gene scoring, no reference)

SC-Type asks a different question: *which canonical marker genes are on/off in each cluster?*

1. **Load functions + database** straight from GitHub via `source()`. The XLSX database lists, per tissue and cell type, the genes you *expect* (positive markers) and *don't expect* (negative markers). We use the "Immune system" tissue. `HGNChelper` fixes outdated gene symbols; `openxlsx` reads the database — that's why those two packages are loaded.
2. **Extract the SCALED matrix.** SC-Type scores on scaled data. The `seurat_v5` check handles the v4/v5 difference in how `scale.data` is stored.
3. **`sctype_score`** produces a cell-type-by-cell score matrix: how strongly each cell matches each candidate type (positive markers raise the score, negative markers lower it).
4. **Aggregate to clusters:** sum scores within each cluster and keep the top-scoring type. We use the known broad labels as the "clusters" here; **in a blind analysis you'd use `seurat_clusters`** — that's the one line you change for real data.
5. **Confidence flag:** clusters whose top score is small relative to their cell count are marked "Low" — a built-in humility check.
6. **Plot** SC-Type calls beside the truth to eyeball agreement.

**Teaching beat:** the default immune database may map several distinct fine states to the *same* SC-Type label (it isn't granular enough). The fix is to **supply your own marker database** in the same XLSX format — which is exactly why the tool accepts a custom DB.

---

## How the three compare
- **Reference mapping** is most accurate *when you have a good, matched reference* — but you need one.
- **SingleR/celldex** is the fastest path when you have no reference of your own; quality depends on atlas match.
- **SC-Type** needs no reference at all, only marker knowledge — great for a first pass or unusual tissues, but only as good as its marker lists.

In practice you often run two of these and reconcile disagreements by hand.

## Check your understanding
1. Why does accuracy fall when you switch from broad to fine labels?
2. SingleR labels per cell while SC-Type labels per cluster — when is each preferable?
3. The default SC-Type database calls two real cell states the same type. What does that tell you, and how do you fix it?
