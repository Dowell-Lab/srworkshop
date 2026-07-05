# Advanced 08 — GSEA over GO Biological Process (gseGO / fgsea)

**Script:** `scripts/advanced/08_gsea_go.R`
**Dataset:** GATA1 (GSE271399), annotated object from `gata1/03`

---

## Why this lesson exists

A differential-expression table gives you a list of genes. But a list of a few hundred genes is hard to interpret and easy to over-read — you cherry-pick the two genes you already know and ignore the rest. **Gene Set Enrichment Analysis (GSEA)** moves the question up a level: instead of asking *"is this one gene significant?"*, it asks *"are the genes belonging to a known pathway collectively shifted up or down?"*

That collective view is far more robust to noise and far easier to turn into a biological story ("erythroid maturation is suppressed", "cell-cycle programs are activated").

---

## The key idea: rank the WHOLE gene list, don't pre-filter

GSEA does not take a filtered list of "significant" genes. It takes **every gene, ranked** from most up-regulated to most down-regulated, and looks at where each pathway's genes fall in that ranking. If a pathway's genes cluster near the top, the pathway is **Activated**; if near the bottom, **Suppressed**.

This is why the DE step in the script runs with **`min.pct = 0` and `logfc.threshold = 0`** — we deliberately keep the full, unfiltered gene list so the ranking is complete.

### The rank statistic

We rank genes by:

```
rank_stat = avg_log2FC * -log10(p_val_adj)
```

This combines **effect size** (fold change) with **confidence** (significance) into one signed number. One subtlety matters: when `p_val_adj == 0` (common with strong effects), `-log10(0)` is infinite. We replace those zeros with `.Machine$double.xmin` (the smallest positive double) so the rank stays finite and the gene keeps its rightful place at the extreme of the ranking.

---

## What the script does, step by step

1. **Load** the annotated GATA1 object and loop over each SingleR cell type.
2. **DE within each cell type**: `FindMarkers(GATA1s vs wtGATA1, min.pct = 0, logfc.threshold = 0)` — the full ranking. (`presto` makes this fast.)
3. **Rank** with `avg_log2FC * -log10(p_val_adj)`, handling `p_val_adj == 0`.
4. **Map SYMBOL → ENTREZID** with `bitr` / `org.Hs.eg.db`, because `gseGO` works in Entrez space.
5. **Run** `gseGO(ont = "BP", by = "fgsea", minGSSize = 15, maxGSSize = 500, eps = 1e-20, pvalueCutoff = 1)` — note `pvalueCutoff = 1` returns **all** terms so *you* filter downstream, and `by = "fgsea"` uses the fast fgsea engine.
6. **Translate** the Entrez IDs in each pathway's `core_enrichment` back to gene symbols so the leading-edge genes are human-readable.
7. **Visualize** — a two-panel NES dot plot (Activated left, Suppressed right; dot size = leading-edge gene count, color = adjusted p-value) plus a volcano of the underlying DE. `pretty_label()` strips the `GOBP_` prefixes for clean plotting.
8. **Write** per-cell-type DE and GSEA CSVs plus a combined table to `OUT_DIR/advanced_gsea_construct/`.

---

## Reading the output

- **NES (Normalized Enrichment Score)** — the signed strength of enrichment. Positive = activated in GATA1s, negative = suppressed. The sign is relative to your contrast direction (`ident.1` = GATA1s).
- **p.adjust** — significance after multiple-testing correction. Filter to `p.adjust < 0.05` before believing a term (the script plots the significant ones).
- **core_enrichment (leading edge)** — the genes actually driving each pathway's signal. This is where the biology lives; read these, not just the pathway name.
- **Two-panel dot plot** — the fastest read: which programs are turned on vs off in each cell type when GATA1s replaces wtGATA1.

---

## Why these parameter choices

- **`by = "fgsea"`** — the fast, permutation-free-ish fgsea algorithm; scales to genome-wide rankings.
- **`minGSSize = 15`, `maxGSSize = 500`** — ignore gene sets that are too tiny (unstable) or too huge (uninformative, e.g. "metabolic process").
- **`eps = 1e-20`** — lets fgsea report very small p-values instead of flooring them, so strongly enriched pathways rank correctly.
- **`pvalueCutoff = 1`** — return everything; filtering at analysis time (not inside `gseGO`) keeps borderline terms visible and lets you re-threshold without re-running.

---

## GSEA vs over-representation (ORA)

A common alternative is **over-representation analysis** (`enrichGO`): take your *significant* genes and ask which pathways are over-represented. GSEA is usually preferable here because it (a) uses the whole ranking so it catches coordinated *small* shifts that never pass a per-gene cutoff, and (b) removes the arbitrary "what counts as significant?" threshold from the pathway step. Use ORA when you genuinely only have a gene list and no ranking.

---

## Common pitfalls
- **Pre-filtering genes before GSEA.** GSEA needs the full ranking — a filtered list breaks the enrichment score. Keep `min.pct = 0`, `logfc.threshold = 0`.
- **Forgetting `p_val_adj == 0`.** `-log10(0)` is infinite and corrupts the ranking; replace with `.Machine$double.xmin`.
- **Mapping problems.** `gseGO` needs Entrez IDs; unmapped or duplicated symbols silently drop genes. Deduplicate after `bitr`.
- **Reading pathway names without the leading edge.** The pathway label is a headline; `core_enrichment` is the story. Always check which genes drove it.
- **Trusting unadjusted or `pvalueCutoff = 1` output as final.** That output is intentionally unfiltered — apply `p.adjust < 0.05` yourself.
- **Comparing NES sign without tracking contrast direction.** Positive NES means "up in `ident.1`"; if you flip the contrast, the signs flip.

## Check your understanding
1. Why does GSEA use the full, unfiltered, ranked gene list instead of just the significant genes?
2. Write the rank statistic used here and explain why `p_val_adj == 0` must be handled specially.
3. Why must gene symbols be converted to Entrez IDs before `gseGO`, and what can go wrong in that step?
4. What do a positive vs negative NES mean for the GATA1s-vs-wtGATA1 contrast, and what determines the sign?
5. When would over-representation analysis (`enrichGO`) be the more appropriate choice than GSEA?
