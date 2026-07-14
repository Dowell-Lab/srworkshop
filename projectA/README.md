# Single-Cell RNA-seq Teaching Course

Authors: Mary Allen and Perplexity AI, based on scripts from Mary Allen, Chris Ozeroff, Jesse Kurland, Georgia Barone, Hope Townsend, Tyler Amos, Brianna Fernandez, Alyx Gray


A hands-on walkthrough of a complete scRNA-seq analysis, from raw counts to
cell-cell communication and developmental trajectories. Each pipeline **step**
has a runnable R script and a companion teaching `.md` lesson explaining what
every block of code does and *why*.

**Audience:** comfortable with R; new to single-cell.

---

## Mapping and counting

In single cell mapping and counting are done in a single step on the super computer. Our mapping example is not on the data we use for the rest of the single cell project. The data for this project is big and expensive so we are not going to map it. We are just showing you how to map using cellranger. 

## Two tracks, two datasets

The course is organized into **two parallel tracks**. You run the GATA1 track
end to end as the primary worked example; the fetal bone marrow track is a
reference you can dip into for a few techniques the GATA1 data doesn't exercise.

### Track A — GATA1 (PRIMARY) · GSE271399

The main dataset you run yourself. An **isogenic iPSC differentiation time
course** of hematopoietic progenitors with a clean 2 × 2 × 3 design:
**genotype** (Euploid / T21) × **construct** (wtGATA1 / GATA1s) × **day**
(D7 / D9 / D11) = 12 conditions. GATA1s is the truncated isoform central to
Down-syndrome-associated leukemia.

- Biology: Takasaki et al., *Single-cell transcriptomics reveal individual and
  cooperative effects of trisomy 21 and GATA1s on hematopoiesis*,
  [*Stem Cell Reports* 2025](https://pubmed.ncbi.nlm.nih.gov/40680731/) ·
- If you want to download the data yourself you can do that here:

  [GEO GSE271399](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271399)
- Why primary: clean, well-designed, and a real **time course**, so pseudotime
  can be validated against sampling day and metadata can be derived cleanly from
  sample names.

### Track B — Fetal bone marrow (REFERENCE) · T21 / D21

Two fetal bone marrow samples used to demonstrate techniques the GATA1 data
doesn't need — and, deliberately, to show what **messy metadata** looks like
(one fetus is labeled with two impossible ages; see the pitfalls note).

| Label | Meaning |
|---|---|
| **T21** | Trisomy 21 (Down syndrome) |
| **D21** | Disomic (typical control) |

Biology source: [fetal bone marrow cell atlas](https://developmental.cellatlas.io/fetal-bone-marrow).

---

## You can run this local or on a super computer. 

If you run it local make sure you have installed the packages in Local_installs.R.
You also need to download the data for the anaysis. You can do this by clicking here. 
Finnaly, you need to change the paths to your data paths. That infomration is in 00_paths_and_setup.R Every single script is going to get the data paths for writing and reading from 00_paths_and_setup.R. So make sure it works before starting. 

## If you work on the super computer. 
We are sharing--- so it might be slow and it might crash. 
Data lives read-only on the cluster — no downloads.
https://ec2-13-56-196-160.us-west-1.compute.amazonaws.com/s/57ea13c286bd33c286bd3/workspaces/

All data is **read-only on the shared AWS/Fiji cluster**. Nobody downloads
anything. Paths are centralized in `scripts/00_paths_and_setup.R`; if the
workshop directory ever moves, that's the only file to edit. Every script writes
its intermediates to a per-user, writable `OUT_DIR` (defaults to
`~/scRNAseq_course_output`), never to the share.

| Track | Stage | Location on the share |
|---|---|---|
| GATA1 | Raw matrices (all lessons) | `…/cookingShow/day7a/iPCStoblood/GSE271399_<sample>_{matrix.mtx,barcodes.tsv,features.tsv}.gz` |
| FBM | Raw 10x matrices (QC) | `…/cookingShow/day7a/total_data/{T21BM_male04, D21_male35}` |
| FBM | Pre-labeled Seurat objects | `…/cookingShow/day9/cellchat-prelabeled-seurat-objs/{t21,d21}_norm_seurat_obj.RData` |
| FBM | Pre-made CellChat objects | `…/cookingShow/day9/cellchatobjs/{t21,d21}-cellchat-obj.rds` |

---

## Course map

### Track A — GATA1 (run these in order)

| # | Lesson | Script | What you learn |
|---|---|---|---|
| 1 | [Load, QC & metadata-from-names](lessons/gata1/01_load_qc_metadata.md) | `scripts/gata1/01_load_qc_metadata.R` | ReadMtx per sample, merge, QC metrics & filtering, **parsing metadata from sample names**, stress/apoptosis module scores |
| 2 | [Clustering, UMAP & composition](lessons/gata1/02_cluster_umap_composition.md) | `scripts/gata1/02_cluster_umap_composition.R` | Normalize/scale/PCA, clustering, UMAP, **composition heatmaps with hierarchical clustering** |
| 3 | [Multi-reference annotation](lessons/gata1/03_multiref_singler_annotation.md) | `scripts/gata1/03_multiref_singler_annotation.R` | JoinLayers, **SingleR across several celldex atlases**, confidence pruning, marker validation |
| 4 | [Cell-cell communication](lessons/gata1/04_cellchat.md) | `scripts/gata1/04_cellchat.R` | CellChat built **from scratch**: object, network inference, pathway analysis |
| 5 | [Pseudotime / trajectories](lessons/gata1/05_pseudotime.md) | `scripts/gata1/05_pseudotime.R` | monocle3: lineage subset, learn_graph, headless root picking, **validating pseudotime against sampling day** |

### Track B — Fetal bone marrow (reference)

| # | Lesson | Script | What you learn |
|---|---|---|---|
| 1 | [QC, Filtering & Integration](lessons/fbm/01_qc_filtering_integration.md) | `scripts/fbm/01_qc_filtering_integration.R` | Read10X, QC, SoupX, RPCA integration, clustering, UMAP |
| 2 | [Marker Genes & Doublets](lessons/fbm/02_markers_doublets.md) | `scripts/fbm/02_markers_doublets.R` | FindAllMarkers/FindMarkers, condition DE, DoubletFinder |
| 3 | [Cell Type Annotation](lessons/fbm/03_cell_type_annotation.md) | `scripts/fbm/03_cell_type_annotation.R` | Seurat reference mapping, SingleR + celldex, SC-Type marker scoring |
| 4 | [Cell-Cell Communication](lessons/fbm/04_cellchat.md) | `scripts/fbm/04_cellchat.R` | CellChat from **pre-made** objects, pathway analysis, condition comparison |
| 5 | [Pseudotime / Trajectories](lessons/fbm/05_pseudotime.md) | `scripts/fbm/05_pseudotime.R` | monocle3: cds setup, learn_graph, order_cells, trajectory-variable genes |

### Cross-cutting teaching note

- [**Metadata pitfalls**](lessons/metadata_pitfalls.md) — contrasts the clean
  "parse from the sample name" approach (GATA1) with the hand-typed sheet that
  put two impossible ages on one FBM fetus. Read after GATA1 Lesson 1.

---

## How to run

1. `cd scripts/gata1/` (or `scripts/fbm/`).
2. Open the script for the lesson you're on; read the matching `lessons/**/…md`
   alongside it.
3. Each script begins with `source("../00_paths_and_setup.R")` to get the
   cluster paths and a writable `OUT_DIR`.
4. Within each track, lessons chain through saved `.rds` objects in `OUT_DIR`, so
   run them in order. The GATA1 track flows 1 → 2 → 3, and lessons 4 and 5 both
   build on the annotated object from lesson 3.

---

## Required packages

| Package | Used in | Purpose |
|---|---|---|
| Seurat | all | core single-cell object & analysis |
| Matrix | GATA1 1 | sparse matrix I/O for ReadMtx |
| tidyverse (dplyr, tidyr, ggplot2) | most | data wrangling & plotting |
| patchwork | several | combining plots |
| SoupX | FBM 1 | ambient RNA correction |
| DoubletFinder | FBM 2 | doublet detection |
| SingleR | GATA1 3, FBM 3 | reference-based annotation |
| celldex | GATA1 3, FBM 3 | built-in reference atlases |
| HGNChelper | FBM 3 | gene-symbol correction (SC-Type) |
| openxlsx | FBM 3 | read SC-Type marker database |
| CellChat | GATA1 4, FBM 4 | cell-cell communication |
| ComplexHeatmap | GATA1 4 | CellChat signaling-role heatmaps |
| monocle3 | GATA1 5, FBM 5 | trajectory / pseudotime |
| R.utils | GATA1 5, FBM 5 | utility functions for monocle3 workflow |
| SeuratWrappers | GATA1 5, FBM 5 | Seurat ↔ monocle3 conversion |

External tools referenced: **CellRanger 7.2.0** produced the raw matrices.

---

## Notes for instructors

- **GATA1 is the spine of the course.** Students run all five GATA1 lessons on
  data they process themselves — including building CellChat from scratch and
  validating pseudotime against the D7/D9/D11 clock.
- **The FBM track adds contrast and pre-computed shortcuts.** Its labeled objects
  carry true cell-type labels at two resolutions (`broad_extfig7A_cell.labels`,
  `cell.labels`) that FBM Lesson 3 uses to *measure* annotation accuracy, and its
  CellChat objects are pre-made so class time can focus on interpretation.
- **Metadata pitfalls is the unifying lesson.** Use the FBM two-ages error to
  motivate the GATA1 parse-from-names discipline.
- Slow steps (CellChat `computeCommunProb`, monocle3 `order_cells`) are clearly
  marked. The GATA1 pseudotime script includes a **headless root picker** so it
  runs without an interactive display.
- `load_one()` (in `00_paths_and_setup.R`) grabs whatever object name a labeled
  `.Rdata` file holds, so the FBM scripts are robust to naming differences.
- Every script writes intermediates to a per-user `OUT_DIR`, never to the
  read-only share.
