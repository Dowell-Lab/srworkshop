# GATA1 Lesson 03 — Multi-reference SingleR annotation and marker checks

**Script:** `scripts/gata1/03_multiref_singler_annotation.R`
**Input:** `OUT_DIR/gata1_combined_clustered.rds` (from Lesson 02)

Clusters from Lesson 02 are just numbered groups. This lesson gives them
biological names automatically with **SingleR**, shows why you might consult
*several* reference atlases, and then verifies the labels against marker genes
you trust. Automated labels are a hypothesis — markers are the check.

## Step 1 — One matrix for SingleR (`JoinLayers`)

After `merge()`, the RNA assay stores a separate `data` layer per sample
(`data.EuploidGATA1sD7`, `data.T21wtGATA1D9`, …). SingleR needs **one** matrix
of all cells, so we call `JoinLayers()` to stitch them into a unified `data`
layer, then pull it out with `LayerData(..., layer = "data")`. We feed the
log-normalized values, which is what SingleR's correlation scoring expects.

## Step 2 — Why multiple references

`SingleR` labels each cell by correlating its profile against a reference of
known cell types ([Aran et al., *Nat Immunol* 2019](https://www.nature.com/articles/s41590-018-0276-y)).
The label you get depends on the reference you pick. `celldex` ships several,
and for iPSC-derived blood each has trade-offs:

| Reference | Strength for this dataset |
|---|---|
| `HumanPrimaryCellAtlasData()` | Broad human primary cells — good first pass (our default) |
| `BlueprintEncodeData()` | Broad stroma + immune |
| `MonacoImmuneData()` | Fine-grained immune subsets |
| `NovershternHematopoieticData()` | Classic hematopoietic hierarchy — blood-specific |
| `DatabaseImmuneCellExpressionData()` | DICE immune populations |

`searchReferences()` lists more. The teaching point: **run a broad reference,
then a tissue-specific one, and compare.** Agreement between independent
references is far stronger evidence than one label alone. The script annotates
with the broad atlas and leaves a commented-out block to add
`NovershternHematopoieticData()` and cross-tabulate the two label sets.

References load from the `celldex` cache — they are **not** part of the
read-only data share.

## Step 3 — Confidence and the long tail

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

We write `SingleR_label`, `SingleR_labels_other`, `SingleR_pruned`, and
`SingleR_delta` back into the Seurat metadata and plot all three label versions
on the UMAP.

## Step 4 — Validate with canonical markers

Never trust an automated label you haven't sanity-checked. `FeaturePlot`s
overlay known genes on the UMAP:

- **`HBG2`, `HBG1`** — fetal hemoglobin → erythroid cells
- **`GATA1`** — the master erythroid/megakaryocyte TF, the gene of interest
- **`GFI1B`** — erythroid/megakaryocyte
- **`TFRC` (CD71)** — erythroid progenitors
- **`CD34`** — hematopoietic stem/progenitor cells
- **`PTPRC` (CD45)** — pan-leukocyte
- **`ITGA4`, `VAMP8`, `DIAPH3`** — additional markers used in the paper

If the cluster SingleR called "erythroid" lights up for `HBG1`/`GATA1`/`TFRC`,
the label holds. If markers and labels disagree, believe the markers and
investigate. The script also saves a side-by-side of labels next to `HBG1`.

The fully annotated object is saved as `gata1_combined_annotated.rds`.

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
