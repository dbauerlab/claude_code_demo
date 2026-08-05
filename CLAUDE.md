# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is a **data-only repository** — it contains an RNA-seq gene-expression dataset and no analysis code, build system, or tests. Any analysis scripts, environment, or pipeline must be created from scratch; do not assume tooling exists.

The dataset is the well-known **"airway" experiment**: human airway smooth muscle cell lines treated with the glucocorticoid **dexamethasone**, used for differential gene expression analysis. Counts were produced by RSEM (as merged by the nf-core/rnaseq pipeline).

## Data files (`data/`)

- **`experiment_table.csv`** — sample sheet. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine` (one of 4: `N61311`, `N052611`, `N080611`, `N061011`), `treatment` (`Untreated` or `Dexamethasone`). 8 samples total = **4 cell lines × 2 conditions**, a paired design.

- **`rsem.merged.gene_counts.tsv`** — tab-separated gene-level counts. Columns: `gene_id` (Ensembl `ENSG` ID), `transcript_id(s)` (comma-separated Ensembl `ENST` IDs), then one column **per sample named by its SRR accession**. ~58,735 gene rows.

## Key facts for working with the data

- **Counts are RSEM estimated counts and are non-integer floats** (e.g. `284.50`). Tools expecting raw integer counts (DESeq2, edgeR) require rounding first (`round()`); do not assume integers.
- **Join key between the two files**: the sample-column headers in the counts TSV match the `ID` column in `experiment_table.csv`. Column order in the TSV is not guaranteed to match row order in the sample sheet — always align by ID, not position.
- **Suggested DE model**: `~ cellLine + treatment` (control for cell line as a paired/blocking factor; test the treatment effect). Set `Untreated` as the reference level.
- The `transcript_id(s)` column must be dropped before building a numeric count matrix.

## Notes

- There is currently a single commit and no README, license, or dependency manifest. When adding analysis, consider recording the environment (e.g. an R `DESeq2` or Python `pydeseq2` setup) so results are reproducible.
