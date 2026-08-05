# CLAUDE.md

This file gives guidance to Claude Code (and other AI assistants) when working in this repository.

## What this project is

`claude_code_demo` is a small bioinformatics repository holding an **RNA-seq gene expression dataset** used for demonstration and analysis. It currently contains data only — there is no analysis code, pipeline, or environment definition checked in yet.

The data is a classic **glucocorticoid (dexamethasone) treatment experiment in airway smooth muscle cells**: four human cell lines, each sequenced under two conditions (untreated vs. dexamethasone-treated), giving 8 samples. It matches the well-known "airway" RNA-seq dataset (Himes et al., 2014), with sample IDs from SRA (`SRR10395xx`).

## Repository layout

```
data/
  experiment_table.csv            # Sample metadata (the design table)
  rsem.merged.gene_counts.tsv     # RSEM merged gene-level counts matrix
```

- **`data/experiment_table.csv`** — 8 samples, columns: `ID`, `cellLine`, `treatment`.
  - `ID`: SRA run accession (e.g. `SRR1039508`)
  - `cellLine`: donor cell line (`N61311`, `N052611`, `N080611`, `N061011`)
  - `treatment`: `Untreated` or `Dexamethasone`
  - The design is paired: each of the 4 cell lines appears once untreated and once treated.

- **`data/rsem.merged.gene_counts.tsv`** — gene-by-sample count matrix produced by RSEM (as merged by the nf-core/rnaseq pipeline).
  - ~58,735 genes (rows) × 8 samples.
  - Columns: `gene_id` (Ensembl gene ID, e.g. `ENSG00000000003`), `transcript_id(s)` (comma-separated Ensembl transcript IDs), then one column per sample named by its `SRR...` accession.
  - Values are estimated (non-integer) expected counts; import with `read.table(..., sep="\t", header=TRUE)` / `pandas.read_csv(sep="\t")`.

The `ID` values in `experiment_table.csv` correspond exactly to the sample column headers in the counts matrix — join on these to attach experimental conditions to expression data.

## Typical use

A natural analysis is **differential expression: Dexamethasone vs. Untreated**, controlling for cell line (paired design). Standard tooling:

- **R / Bioconductor**: DESeq2 or edgeR. With DESeq2, round the RSEM counts to integers, build the design as `~ cellLine + treatment`, and set `treatment` reference level to `Untreated`.
- **Python**: pandas for loading, `pydeseq2` for the DE model.

Genes with all-zero counts (e.g. `ENSG00000000005`) should be filtered before modeling.

## Repository info

- Remote: `https://github.com/dbauerlab/claude_code_demo`
- Default branch: `main`
- Lab: Bauer Lab (`dbauerlab`)

## Notes for assistants

- The counts TSV is ~7 MB; prefer streaming / column-aware loading over dumping the whole file into context.
- There is no build system, test suite, or dependency manifest here yet. If you add analysis code, also add an environment file (`environment.yml`, `renv`, or `requirements.txt`) so the analysis is reproducible.
- Keep raw data in `data/` unmodified; write derived outputs to a separate directory (e.g. `results/`).
