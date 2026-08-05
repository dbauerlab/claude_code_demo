# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository overview

This is a data repository for an RNA-seq differential expression analysis. It contains no analysis scripts yet — only the raw input data files.

## Data

All data lives in `data/`:

- **`experiment_table.csv`** — sample metadata. Columns: `ID` (SRR accession), `cellLine` (4 human airway smooth muscle cell lines: N61311, N052611, N080611, N061011), `treatment` (Untreated or Dexamethasone). 8 samples total: one untreated/treated pair per cell line.

- **`rsem.merged.gene_counts.tsv`** — RSEM gene-level count matrix. Rows are genes (~58,735 Ensembl gene IDs), columns are `gene_id`, `transcript_id(s)`, then one column per SRR sample. Values are raw counts (floats from RSEM). Sample column order matches the `ID` column in `experiment_table.csv`.

## Dataset context

This is the classic [airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html) dataset (Himes et al. 2014, PubMed PMID: 24926665) — commonly used in Bioconductor DESeq2/edgeR tutorials. The biological question is the effect of dexamethasone (a glucocorticoid) on gene expression in human airway smooth muscle cells.
