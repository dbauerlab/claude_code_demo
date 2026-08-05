# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

This is a **data-only repository** — there is no application code, build system, or test
suite. It contains a bulk RNA-seq gene-count dataset intended for downstream differential
expression analysis. Any code (analysis scripts, notebooks) does not yet exist and would be
added by the user; there are no project-defined commands to build, lint, or test.

The dataset is the well-known **"airway" experiment**: primary human airway smooth muscle
cell lines treated with dexamethasone (a glucocorticoid) versus untreated controls.

## Data layout

Both files live in `data/` and are joined on the sample IDs (SRR accessions).

- `data/experiment_table.csv` — the sample sheet / experimental design. 8 samples, columns:
  - `ID` — SRA run accession (e.g. `SRR1039508`), matches the count-matrix column headers
  - `cellLine` — donor cell line (4 lines: N61311, N052611, N080611, N061011)
  - `treatment` — `Untreated` or `Dexamethasone`
  - The design is **paired**: each of the 4 cell lines contributes one untreated and one
    treated sample. Account for `cellLine` as a blocking factor when modelling
    treatment effects (see Analysis conventions).

- `data/rsem.merged.gene_counts.tsv` — RSEM-generated merged gene-level count matrix
  (~58,735 genes × 8 samples), tab-separated. Columns:
  - `gene_id` — Ensembl gene ID (`ENSG...`)
  - `transcript_id(s)` — comma-separated Ensembl transcript IDs collapsed into the gene
  - one column per sample, headed by the SRA accession, holding **estimated counts**
    (RSEM expected counts, hence non-integer / fractional values)

## Working with this data

- Counts are RSEM **expected counts** (fractional). DESeq2 (and edgeR) require integer
  counts — **round the matrix to integers before importing** (e.g.
  `round(as.matrix(counts))`), unless using a method built for estimated counts
  (e.g. `tximport`).
- Match count-matrix sample columns to `experiment_table.csv` rows via the SRA accession.
- The count matrix is ~6.5 MB; prefer streaming/columnar reads over loading naively when
  only a subset of samples or genes is needed.
- `gene_id` values are Ensembl IDs without version suffixes.

## Analysis conventions

- **Language: R.** Analysis is done in R (DESeq2 for the DE workflow). Assume Bioconductor
  packages unless told otherwise.
- **Model:** `~ cellLine + treatment` — `treatment` is the effect of interest, `cellLine`
  is the paired blocking factor controlling for donor.
- **Reference level:** `Untreated` is the reference. Set it explicitly so log2 fold changes
  read as Dexamethasone-vs-Untreated:
  `colData$treatment <- relevel(factor(colData$treatment), ref = "Untreated")`.
  Do the same for `cellLine` if its ordering matters.
- **Deliverables are Quarto `.qmd` documents.** Write analyses as `.qmd` files (rendered with
  Quarto), not bare `.R` scripts or `.Rmd`, so narrative, code, and figures ship together.
