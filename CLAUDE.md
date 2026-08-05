# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics demo project for RNA-seq differential expression analysis of the *airway* dataset. Apart from this file the repository is **data-only** — the two files in `data/` are the sole other tracked contents, and no analysis code exists yet.

## Language and deliverables

- **All analysis is written in R.** Do not introduce Python for the analysis path.
- **Deliverables are Quarto documents (`.qmd`)**, not bare `.R` scripts. New analysis work belongs in a `.qmd` that renders to a readable report with narrative prose, the code that produced each result, and its figures and tables inline.
- Render with `quarto render <file>.qmd`. There is no dependency manifest, lint config, or test suite yet, so a clean render is the de facto check that an analysis still works — render before considering a change done.
- Core toolchain is Bioconductor: DESeq2 for the DE model, plus the usual tidyverse for data handling and plotting.

## Data

- **`data/experiment_table.csv`** (268 B) — sample metadata, columns `ID,cellLine,treatment`. 8 samples: 4 human airway smooth muscle cell lines (N61311, N052611, N080611, N061011), each contributing one `Untreated` and one `Dexamethasone` sample. `ID` values are NCBI SRA run accessions.
- **`data/rsem.merged.gene_counts.tsv`** (6.5 MB) — gene-level count matrix from RSEM. Columns: `gene_id` (Ensembl, unversioned), `transcript_id(s)` (comma-separated), then one column per sample named by SRA accession.

This is the classic **airway** dataset (Himes et al. 2014, PMID 24926665), the standard worked example in DESeq2 tutorials. Published results are therefore a useful sanity check — dexamethasone should strongly induce *DUSP1*, *KLF15*, *PER1*, and *CRISPLD2*.

## Data characteristics that affect analysis

Measured from the actual files — worth knowing before writing loading code:

- **58,735 genes**, of which **28,344 (48%) are all-zero across every sample.** Pre-filter (e.g. `rowSums(counts) > 0`, or a stricter minimum-count threshold) before DE testing.
- **Library sizes range 14.9 M–30.3 M** assigned reads (smallest SRR1039513, largest SRR1039517) — a ~2× spread, so rely on DESeq2's size-factor normalisation rather than raw counts.
- **No gene symbol column.** Mapping Ensembl IDs to symbols needs external annotation (`org.Hs.eg.db`, biomaRt, or an `annotables`-style table).

## Experimental design and analysis conventions

These are fixed project conventions. Follow them in any new analysis rather than re-deriving a design.

**The design is paired by donor cell line.** Each of the 4 cell lines provides its own untreated control, so the cell line is a blocking factor and must stay in the model:

```r
design = ~ cellLine + treatment
```

`treatment` is the variable of interest; `cellLine` absorbs donor-to-donor baseline differences. Dropping `cellLine` discards the pairing and inflates within-group variance, costing real power on only 8 samples.

**`Untreated` is the reference level, and must be set explicitly.** R orders factor levels alphabetically, which would make `Dexamethasone` the reference and silently invert the sign of every log2 fold change:

```r
coldata$treatment <- relevel(factor(coldata$treatment), ref = "Untreated")
coldata$cellLine  <- factor(coldata$cellLine)
```

With this reference, a positive LFC means dexamethasone *induces* the gene.

**RSEM counts are non-integer and must be rounded.** RSEM apportions multi-mapping reads probabilistically, so 11,117 genes carry fractional values. DESeq2 requires integers and will error otherwise:

```r
counts <- round(as.matrix(counts))
```

Use `round()`, not `as.integer()`/`floor()`/`trunc()` — those truncate downward and bias low-count genes toward zero.

**Sample order is not guaranteed.** Subset the count columns by the metadata IDs (`counts[, coldata$ID]`) so columns and rows correspond by name, not by incidental file ordering.
