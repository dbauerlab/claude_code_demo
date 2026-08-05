# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository state

This repository currently contains only input data — no analysis code, scripts, notebooks, build tooling, or tests exist yet. There is no README. Treat any request here as starting a bioinformatics analysis from scratch on top of these two files.

## Data

- `data/experiment_table.csv` — sample sheet, 8 samples. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine` (donor identifier, e.g. `N61311`), `treatment` (`Untreated` or `Dexamethasone`). Samples are paired: each of the 4 cell lines has one Untreated and one Dexamethasone run, making cell line a natural blocking/pairing factor for differential expression (paired design by donor).
- `data/rsem.merged.gene_counts.tsv` — gene-level expression matrix from RSEM, merged across samples (the format nf-core/rnaseq's RSEM merge step produces). Columns: `gene_id` (Ensembl gene ID, e.g. `ENSG00000000003`), `transcript_id(s)` (comma-separated Ensembl transcript IDs for that gene), then one column per sample matching the `ID` values in `experiment_table.csv`. ~58,735 genes.
- Values in the count matrix are RSEM **expected counts**, which are non-integer floats (e.g. `340.80`). Downstream tools that require integer counts (DESeq2, edgeR) need these rounded first.
- This is the well-known "airway" dexamethasone RNA-seq dataset (Himes et al.) — airway smooth muscle cells from 4 donors, treated with dexamethasone or left untreated.

## Working in this repo

- Sample identity is joined between the two files via `ID` in `experiment_table.csv` matching the sample column headers in `rsem.merged.gene_counts.tsv`.
- Since there's no existing code structure or dependency manifest, when adding analysis code, ask the user which language/tool stack they want (e.g. R/DESeq2, Python/pydeseq2) rather than assuming.
