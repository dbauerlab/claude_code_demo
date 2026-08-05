# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is a **data repository**, not a software project. It currently contains no source code, build system, tests, or dependencies — only bulk RNA-seq quantification data and its experimental design. Any analysis code (differential expression, QC, plotting) would need to be added.

## Contents

- `data/experiment_table.csv` — sample sheet / experimental design. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine`, `treatment`. 8 samples across 4 cell lines (`N61311`, `N052611`, `N080611`, `N061011`), each with a paired `Untreated` / `Dexamethasone` condition.
- `data/rsem.merged.gene_counts.tsv` — RSEM merged gene-level counts, tab-separated, ~58.7k genes (rows) × 8 samples. Columns: `gene_id` (Ensembl gene ID, e.g. `ENSG00000000003`), `transcript_id(s)` (comma-separated Ensembl transcript IDs mapping to the gene), then one count column per sample named by SRA accession. Counts are RSEM estimates, so values are floats, not integers.

## Key facts for analysis

- The sample columns in `rsem.merged.gene_counts.tsv` join to `experiment_table.csv` on the SRA accession (`ID` column ↔ count column headers). Confirm the 8 accessions match before analysis.
- This is the classic **airway** dataset: dexamethasone treatment of human airway smooth muscle cell lines.

## Experimental design & analysis conventions

- **Language:** All analysis is done in **R**.
- **Deliverables:** Analyses are authored as **Quarto (`.qmd`) documents**, not plain scripts.
- **Design:** Paired 4×2 (cell line × treatment). Use the model formula **`~ cellLine + treatment`** — `cellLine` is a blocking factor controlling for the paired structure, and `treatment` is the effect of interest.
- **Reference level:** `Untreated` is the reference (baseline) level for `treatment`; relevel the factor explicitly (e.g. `relevel(treatment, ref = "Untreated")`) so log-fold-changes read as Dexamethasone-vs-Untreated.
- **Counts are non-integer:** RSEM counts are floats. DESeq2 requires integer counts, so **round before constructing the `DESeqDataSet`** (e.g. `round(counts)`). RSEM-aware import via `tximport` is the more correct alternative when transcript-level files are available.
