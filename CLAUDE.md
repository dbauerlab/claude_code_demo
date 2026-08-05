# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bulk RNA-seq differential expression demo. The data is a dexamethasone treatment experiment in human airway smooth muscle cells (the [Himes et al. 2014 / `airway` dataset](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4057696/)). The analysis goal is differential gene expression, Dexamethasone-treated vs Untreated.

## Repository State

Only the two input data files are tracked; there are no analysis scripts, no package manifest (no `renv.lock`, `DESCRIPTION`, `environment.yml`, `requirements.txt`), and no CI. **There are therefore no build, lint, or test commands to run** — don't go looking for them. The language and environment for the analysis are still an open choice; ask before assuming R/Bioconductor vs Python.

## Data

**`data/experiment_table.csv`** — sample metadata, 8 rows:

| column | meaning |
| --- | --- |
| `ID` | SRA run accession — the sample identifier used everywhere |
| `cellLine` | donor: `N61311`, `N052611`, `N080611`, `N061011` |
| `treatment` | `Untreated` or `Dexamethasone` |

Fully balanced paired design: each of the 4 donor cell lines contributes exactly one untreated and one treated sample. There are no unpaired samples and no other covariates.

**`data/rsem.merged.gene_counts.tsv`** — RSEM gene-level counts, 58,735 gene rows × 8 sample columns (~6.8 MB):

- `gene_id` — Ensembl gene ID; all rows match `^ENSG`, none carry a `.version` suffix, and there are no duplicates, so it is safe to use directly as a row key.
- `transcript_id(s)` — **column 2 is not a sample.** It is a comma-separated transcript list. Code that treats every column after the first as a count column will silently produce a garbage sample. Drop it explicitly.
- Columns 3–10 are the 8 samples, and they happen to appear in the same order as the metadata rows. Match by name anyway rather than relying on position.

## Key Analysis Considerations

- **Block on cell line.** The paired design means `cellLine` belongs in the model as a blocking factor (`~ cellLine + treatment` in DESeq2). Omitting it leaves donor-to-donor variation in the residuals and costs real power.
- **Round before count-based tools.** RSEM apportions multi-mapping reads probabilistically, so values are fractional — 11,117 rows contain at least one non-integer. DESeq2 and edgeR require integers.
- **Pre-filter aggressively.** 28,344 of the 58,735 genes (~48%) are zero in all 8 samples. A total-count ≥ 10 filter leaves ~20,783 genes.
- **Library sizes vary about 2×** (14.9M–30.3M summed counts, smallest `SRR1039513`, largest `SRR1039517`), so work from normalized values for any visual or distance-based comparison.
- **Gene symbols need mapping.** IDs are bare Ensembl; map via `org.Hs.eg.db` or `biomaRt`. Expect unmapped IDs and a non-unique symbol mapping — decide how to handle collisions rather than dropping rows silently.
