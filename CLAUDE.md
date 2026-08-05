# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **bioinformatics RNA-seq dataset**, currently data-only — there is no analysis code, build system, or test suite yet. It holds the inputs for a differential gene expression analysis and is expected to grow analysis scripts (likely R/DESeq2 or Python) over time.

The data is the well-known **airway** dataset (Himes et al. 2014): human airway smooth muscle cells, comparing dexamethasone (a glucocorticoid) treatment against untreated controls, across four cell-line donors. This is a paired design — each donor contributes one treated and one untreated sample.

## Data layout

Everything lives in `data/`:

- **`experiment_table.csv`** — the sample metadata / experimental design. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine` (donor: `N61311`, `N052611`, `N080611`, `N061011`), `treatment` (`Untreated` or `Dexamethasone`). 8 samples = 4 donors × 2 conditions.
- **`rsem.merged.gene_counts.tsv`** — the RSEM-merged gene-level count matrix. Tab-separated. First two columns are `gene_id` (Ensembl `ENSG…`) and `transcript_id(s)` (comma-separated `ENST…`); the remaining 8 columns are one per sample, keyed by the same `ID` accession as `experiment_table.csv`. ~58,735 genes.

## Facts worth knowing before analysis

- **The two files join on the sample ID.** The count matrix's sample column headers (`SRR…`) are exactly the `ID` values in `experiment_table.csv` — this is the key that links counts to their treatment/donor labels.
- **Counts are RSEM expected counts, not raw integers.** Values are frequently fractional (e.g. `284.50`) because RSEM distributes multi-mapping reads probabilistically. Tools expecting integer counts (e.g. DESeq2) need these rounded first.
- **The design is paired on `cellLine`.** A correct DE model should control for donor — i.e. `~ cellLine + treatment`, not `~ treatment` alone.
