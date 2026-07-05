# Advanced 06 — Harmony integration + LISI/cLISI diagnostics

**Script:** `scripts/advanced/06_harmony_integration_lisi.R`
**Dataset:** fetal bone marrow (T21 + D21 pre-labeled objects)

---

## Why this lesson exists

The basic FBM track integrated samples with **RPCA** and then eyeballed the UMAP to decide whether it "looked integrated." That is fine for a first pass, but it is subjective and it doesn't tell you whether you *over*-corrected — that is, whether you accidentally erased real biology while chasing batch mixing.

This lesson teaches two things the course did not yet cover:

1. **Harmony** — a fast, widely used alternative integration method that corrects batch effects directly in PCA space.
2. **LISI** — a *quantitative* score that tells you whether integration actually worked, so the decision is measured rather than a vibe.

---

## The core tension in every integration

Every integration is a tug-of-war between two goals:

- **Mix the batches** — cells that differ only because they came from different samples/donors should end up sitting together.
- **Preserve the cell types** — cells that differ because they are genuinely different cell types should stay apart.

Push too hard on mixing and you smear distinct cell types into each other. Push too little and batch effects dominate your clusters. You need a number for each side of that tug-of-war.

---

## What LISI measures

**LISI = Local Inverse Simpson's Index.** For each cell, it looks at the neighborhood in the embedding and asks *how many different categories are represented there*, as an effective count.

We compute it two ways (this is the key idea):

- **iLISI** — LISI over the **batch/sample** variable. It answers *"in a cell's neighborhood, how many different samples are present?"* **Higher is better** — high iLISI means samples are well mixed.
- **cLISI** — LISI over the **cell-type label** variable. It answers *"in a cell's neighborhood, how many different cell types are present?"* **Lower is better** — low cLISI (near 1) means each neighborhood is a single, clean cell type.

So a good integration gives you **high iLISI and low cLISI at the same time**. That is exactly the tradeoff plot the script produces.

---

## Harmony's `theta` knob

Harmony's most important parameter is **`theta`**, the diversity-penalty. Bigger `theta` pushes harder to mix batches.

- `theta = 0` → effectively no batch correction (baseline).
- larger `theta` → more aggressive mixing → iLISI goes up, but cLISI can also creep up (cell types start blending).

Because there is no universally correct value, the script **sweeps** `thetas = c(0, 1, 2, 3, 5)`, computes iLISI and cLISI at each, and plots the tradeoff. You choose the `theta` that maximizes mixing *without* letting cell-type LISI climb — the classic "elbow" judgment call, now backed by numbers.

---

## What the script does, step by step

1. **Build one merged object** from the two labeled FBM objects using `load_one()`, stamping a `$sample` (batch) and `$condition` (T21/D21) column so we have explicit batch and biology variables.
2. **Standard preprocessing** — `NormalizeData` → `FindVariableFeatures(2000)` → `ScaleData` → `RunPCA(npcs = 20)`.
3. **Baseline** — UMAP + clustering + LISI on plain PCA (no correction).
4. **Per-theta loop** — `RunHarmony(group.by.vars = "sample", theta = t)` for each `theta`, then UMAP/cluster/LISI on the `harmony` reduction.
5. **Compare** — assemble iLISI and cLISI across all thetas into a **tradeoff plot** and per-theta UMAPs colored by sample and by cell type.

All reductions and plots are written to `OUT_DIR/advanced_harmony_integration/`. Inputs stay read-only.

---

## Reading the output

- **Baseline UMAP colored by sample** should show batch-driven separation (that's the problem you're fixing).
- **Harmony UMAP colored by sample** should show the samples interleaved.
- **Harmony UMAP colored by cell type** should still show clean, separated cell-type islands.
- **Tradeoff plot** — pick the `theta` on the good side of the elbow: iLISI as high as you can get before cLISI starts rising.

---

## Common pitfalls
- **Confusing the two LISIs.** iLISI (batch) you want *high*; cLISI (cell type) you want *low*. Mixing these up leads to exactly the wrong `theta` choice.
- **Chasing maximum mixing.** The highest iLISI is not the goal — it usually means you've over-corrected and destroyed cell-type structure (cLISI climbs).
- **Computing LISI on the wrong space.** Compute it on the embedding you actually use downstream (here the UMAP off the Harmony reduction), not on raw PCA.
- **Forgetting to set the seed.** Harmony and UMAP are stochastic; without `harmony_seed` your sweep isn't reproducible.
- **Running Harmony without joined layers / a single normalized object.** Batch variables must be per-cell metadata on one object, not separate objects.
- **Using too few cells per batch.** LISI neighborhoods become unstable when a sample contributes very few cells.

## Check your understanding
1. In one sentence each, what do iLISI and cLISI measure, and which direction is "good" for each?
2. What does increasing Harmony's `theta` do, and what is the risk of setting it too high?
3. Why does the script sweep several `theta` values instead of picking one?
4. Why is it not enough to just look at a UMAP and say "that looks integrated"?
5. If iLISI is high but cLISI is also high, what has probably happened to your data?
