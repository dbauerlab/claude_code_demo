# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **bioinformatics dataset** for RNA-seq differential-expression analysis. There is currently no analysis code — the repo contains only raw input data under `data/`. Any pipeline, scripts, or notebooks are still to be written.

The data is the well-known **"airway"** experiment: human airway smooth muscle cell lines treated with dexamethasone (a glucocorticoid) versus untreated controls.

## Data files

- `data/experiment_table.csv` — the sample sheet. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine` (one of four donor lines: N61311, N052611, N080611, N061011), `treatment` (`Untreated` or `Dexamethasone`). 8 samples total.
- `data/rsem.merged.gene_counts.tsv` — gene-level expression counts, tab-separated, produced by RSEM (characteristic of nf-core/rnaseq output). ~58,735 rows. Columns: `gene_id` (Ensembl gene ID, e.g. `ENSG00000000003`), `transcript_id(s)` (comma-separated Ensembl transcript IDs), then one count column per sample keyed by the SRA `ID`.

## Experimental design (important for analysis)

The design is **paired**: each of the 4 cell lines contributes one untreated and one dexamethasone-treated sample. Differential-expression models must account for `cellLine` as a blocking factor and test the `treatment` effect, using the design formula `~ cellLine + treatment` (DESeq2/edgeR/limma). The `ID` column in the sample sheet is the join key to the count-matrix column headers.

- **Reference level:** `Untreated` is the reference level for `treatment`, so the `treatment` coefficient reports Dexamethasone-vs-Untreated log2 fold changes. Set it explicitly (e.g. `relevel(treatment, ref = "Untreated")` or `factor(..., levels = c("Untreated", "Dexamethasone"))`) — do not rely on alphabetical default.

## Analysis conventions

- **Language:** analysis is done in **R** (Bioconductor stack — DESeq2 and friends).
- **Deliverables:** reproducible reports as **Quarto `.qmd`** files.
- Counts are RSEM "expected counts" and are non-integer; **round to integers** before feeding DESeq2 (which requires integer counts).
- Gene and transcript identifiers are unversioned Ensembl IDs; map to gene symbols via an Ensembl/biomaRt annotation matching the genome build used to generate the counts.
