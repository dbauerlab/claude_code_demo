# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics demo project for RNA-seq differential expression analysis. The repository is **data-only** — the two files in `data/` are the sole tracked contents. There is no analysis code, dependency manifest, README, or CI yet, so there are no build, lint, or test commands to run. Expect to establish those conventions when the first script lands.

## Data

- **`data/experiment_table.csv`** (268 B) — sample metadata, columns `ID,cellLine,treatment`. 8 samples: 4 human airway smooth muscle cell lines (N61311, N052611, N080611, N061011), each contributing one `Untreated` and one `Dexamethasone` sample. `ID` values are NCBI SRA run accessions.
- **`data/rsem.merged.gene_counts.tsv`** (6.5 MB) — gene-level count matrix from RSEM. Columns: `gene_id` (Ensembl, unversioned), `transcript_id(s)` (comma-separated), then one column per sample named by SRA accession.

This is the classic **airway** dataset (Himes et al. 2014, PMID 24926665), the standard worked example in DESeq2 tutorials. Published results are therefore a useful sanity check — dexamethasone should strongly induce *DUSP1*, *KLF15*, *PER1*, and *CRISPLD2*.

## Data characteristics that affect analysis

Measured from the actual files — worth knowing before writing loading code:

- **58,735 genes**, of which **28,344 (48%) are all-zero across every sample.** Pre-filter (e.g. `rowSums(counts) > 0`, or a stricter minimum-count threshold) before DE testing.
- **11,117 genes carry fractional counts** — RSEM distributes multi-mapping reads probabilistically. DESeq2 and edgeR require integers, so `round()` the matrix before constructing a `DESeqDataSet`. Do not silently truncate with `as.integer()`, which floors.
- **Library sizes range 14.9 M–30.3 M** assigned reads (smallest SRR1039513, largest SRR1039517) — a ~2× spread, so rely on the tools' own size-factor normalisation rather than raw counts.
- **No gene symbol column.** Mapping Ensembl IDs to symbols needs external annotation (`org.Hs.eg.db`, biomaRt, or an `annotables`-style table).
- Sample columns in the count matrix appear in the **same order** as the metadata rows, but DESeq2 matches by name and errors on mismatch — pass `counts[, metadata$ID]` explicitly rather than relying on incidental ordering.

## Intended analysis design

The experiment is **paired by donor cell line**, which is the key modelling point: each of the 4 lines provides its own untreated control. Use

```
~ cellLine + treatment
```

with `treatment` as the variable of interest and `Untreated` as the reference level (set it explicitly — R's default factor ordering would make `Dexamethasone` the reference alphabetically). Dropping the `cellLine` term discards the pairing and inflates within-group variance.
