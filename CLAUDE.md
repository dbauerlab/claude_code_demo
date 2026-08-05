# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **bioinformatics RNA-seq project** for differential gene expression analysis. Deliverables are **Quarto documents (`.qmd`)** that render the analysis narrative alongside its code and figures. The primary workflow is [`airway_de.qmd`](airway_de.qmd); one-shot toolchain install is [`scripts/setup.sh`](scripts/setup.sh).

**DE engine:** the differential-expression statistics are computed by the **clawbio `rnaseq-de` skill (PyDESeq2)** — a re-implementation of the DESeq2 method — invoked from the `.qmd`. **R is the reporting layer only** (input prep, ingesting clawbio's results, all visualization). Downstream pathway enrichment tries clawbio's `pathway-enricher` skill and falls back to `gseapy.enrichr` ([`scripts/enrich_gseapy.py`](scripts/enrich_gseapy.py)) — clawbio's Enrichr upload is form-encoded and 400s on large gene lists, so gseapy (correct multipart API) is the working path. Because PyDESeq2 implements the DESeq2 method, every statistical convention below still applies unchanged.

The data is the well-known **airway** dataset (Himes et al. 2014): human airway smooth muscle cells, comparing dexamethasone (a glucocorticoid) treatment against untreated controls, across four cell-line donors. This is a paired design — each donor contributes one treated and one untreated sample.

## Data layout

Everything lives in `data/`:

- **`experiment_table.csv`** — the sample metadata / experimental design. Columns: `ID` (SRA run accession, e.g. `SRR1039508`), `cellLine` (donor: `N61311`, `N052611`, `N080611`, `N061011`), `treatment` (`Untreated` or `Dexamethasone`). 8 samples = 4 donors × 2 conditions.
- **`rsem.merged.gene_counts.tsv`** — the RSEM-merged gene-level count matrix. Tab-separated. First two columns are `gene_id` (Ensembl `ENSG…`) and `transcript_id(s)` (comma-separated `ENST…`); the remaining 8 columns are one per sample, keyed by the same `ID` accession as `experiment_table.csv`. ~58,735 genes.

## Analysis conventions

- **Deliverables are rendered `.qmd`.** Each deliverable is a **Quarto `.qmd`** file (knitr/R engine) that is rendered to produce the report. Prefer this literate-document workflow over loose scripts. The DE step shells out to the clawbio Python skill from within the document.
- **Design formula: `~ cellLine + treatment`.** The design is paired — each of the four donors contributes one treated and one untreated sample — so the model must control for donor (`cellLine`) rather than using `~ treatment` alone.
- **`Untreated` is the reference level.** Encoded via the contrast direction `treatment,Dexamethasone,Untreated`, so the reported log2 fold-change is Dexamethasone vs. Untreated. Don't rely on default alphabetical factor ordering.
- **Round RSEM counts before DE.** Counts are RSEM *expected counts* and frequently fractional (e.g. `284.50`), because RSEM distributes multi-mapping reads probabilistically. PyDESeq2 (like DESeq2) requires integers — round the matrix before running the DE step.
- **Python 3.11–3.12 for clawbio.** The clawbio stack (pydeseq2/numpy/scipy) needs Python 3.11–3.12; the `.qmd` calls it via a `python3.12` interpreter installed by `scripts/setup.sh`.

## Facts worth knowing before analysis

- **The two files join on the sample ID.** The count matrix's sample column headers (`SRR…`) are exactly the `ID` values in `experiment_table.csv` — this is the key that links counts to their treatment/donor labels.
- **Column order isn't guaranteed to match.** Before constructing the `DESeqDataSet`, align the count-matrix columns to the metadata rows by `ID` rather than assuming they're in the same order.
