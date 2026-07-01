# Teaching note — Metadata pitfalls (and how to avoid them)

This note ties the two course tracks together around a single, career-defining
lesson: **your analysis is only as trustworthy as your metadata.** The two
datasets were chosen partly to contrast a clean approach with a messy one.

## The tale of two datasets

| | GATA1 track (GSE271399) | Fetal bone marrow track (T21/D21) |
|---|---|---|
| Metadata source | **Parsed from sample names in code** | Hand-assembled sample sheet |
| Design encoding | `T21wtGATA1D9` → genotype/construct/day | Free-text fields per sample |
| A known problem | none — one source of truth | **one fetus labeled with two different ages** |
| Fixability | rename a file, everything updates | must hunt down and reconcile by hand |

The fetal bone marrow data contains a genuine, instructive error: **a single
fetus is annotated with two different ages** — which is biologically impossible,
since a given fetus has exactly one gestational age at collection. That kind of
contradiction is the tip of an iceberg. If one field is wrong, you cannot trust
that the others are right, and every downstream comparison built on that
grouping is suspect.

## Why this matters more in single-cell than almost anywhere else

In bulk RNA-seq a metadata error affects a handful of samples. In single-cell,
one mislabeled sample silently mis-stamps **tens of thousands of cells**. Those
cells then flow into clustering, composition tests, differential expression,
CellChat groupings, and pseudotime roots. The error doesn't announce itself — it
just quietly biases every result. A p-value computed on a wrong grouping still
looks like a p-value.

## The clean pattern: name is the single source of truth

The GATA1 track encodes the entire experimental design *in the sample name* and
parses it back out with code (see `lessons/gata1/01_load_qc_metadata.md`):

```r
combined$genotype  <- ifelse(grepl("^Euploid", combined$sample), "Euploid", "T21")
combined$construct <- case_when(grepl("wtGATA1", combined$sample) ~ "wtGATA1",
                                grepl("GATA1s",  combined$sample) ~ "GATA1s")
combined$day       <- case_when(grepl("D7$",  combined$sample) ~ "D7",
                                grepl("D9$",  combined$sample) ~ "D9",
                                grepl("D11$", combined$sample) ~ "D11")
stopifnot(!any(is.na(combined$genotype)))   # fail loudly if a name didn't match
```

Why this is robust:

1. **One source of truth.** The name *is* the metadata. There's no second sheet
   to fall out of sync, so a sample physically cannot carry two contradictory
   ages.
2. **Self-documenting.** Anyone reading `T21wtGATA1D9` knows exactly what it is.
3. **Reproducible.** Re-run the code, get identical metadata every time — no
   manual copy-paste step to get wrong.
4. **Fails loudly, not silently.** The `stopifnot()` checks turn a typo in the
   sample list into an immediate, obvious error instead of a silent `NA` that
   propagates.

## The catch: this only works if names are disciplined

Programmatic parsing pushes the responsibility upstream to **naming
conventions**. Adopt these habits when you generate data:

- Encode every experimental variable in the name, in a fixed order and
  vocabulary: `{genotype}{construct}{day}`.
- Use consistent tokens (`D7`, not `d7`, `Day7`, `7d`).
- Avoid characters that break file globs or R parsing (spaces, `/`, `.`).
- Keep a short data dictionary documenting the scheme.
- Anchor your regex (`^Euploid`, `D7$`) so `D7` doesn't accidentally match
  inside another token.

## When you inherit messy metadata (the FBM reality)

You won't always control naming. When you receive a dataset with a hand-typed
sheet:

1. **Validate before you analyze.** Cross-tabulate fields
   (`table(fetus_id, age)`) and look for impossible combinations — like one
   fetus with two ages. Check that every categorical field has the expected
   number of levels.
2. **Reconcile against a source of truth.** Go back to the wet-lab notebook,
   submission record, or GEO/SRA entry. Do not guess.
3. **Document the fix.** Record what was wrong and how you resolved it, in code,
   so the correction is reproducible and auditable.
4. **Re-derive from names where possible.** If the filenames are more reliable
   than the sheet, parse from them and treat the sheet as something to check
   *against*, not to trust blindly.
5. **When in doubt, exclude.** A dropped ambiguous sample is honest; a wrongly
   labeled one silently corrupts your conclusions.

## The one-sentence takeaway

**Make the data structure the metadata, validate every field for impossible
combinations, and fail loudly — because in single-cell, a metadata bug isn't a
typo, it's a wrong result with a confident p-value.**

## Check your understanding

1. Why is a metadata error more damaging in single-cell than in bulk RNA-seq?
2. What single property of the "parse from the name" approach makes it impossible
   for one sample to have two contradictory ages?
3. What does `stopifnot(!any(is.na(...)))` protect you from, and why is failing
   loudly better than a silent `NA`?
4. You inherit a sample sheet where one donor appears with two sexes. List the
   steps you'd take before running any analysis.
5. Programmatic parsing moves responsibility "upstream." Upstream to what, and
   what habits make it safe?
