# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A **data-only starting point** for a bulk RNA-seq differential expression demo (`dbauerlab/claude_code_demo`). At present the repo contains nothing but `data/` — no source code, no build system, no tests, no dependency manifest, and no README. Any analysis code, and the choice of language and tooling for it, is still open; do not assume an existing R or Python project layout.

Each workshop participant works on their own `demo/<name>` branch off `main`. Keep new work on the current `demo/*` branch rather than committing to `main`.

## The dataset

This is the **airway** dataset (Himes et al. — human airway smooth muscle cells, dexamethasone vs. untreated), quantified with RSEM and merged by the nf-core/rnaseq pipeline.

`data/experiment_table.csv` — the sample sheet, and the only place the experimental design lives:

- 8 samples, columns `ID, cellLine, treatment`
- 4 cell lines (`N61311`, `N052611`, `N080611`, `N061011`), each contributing one `Untreated` and one `Dexamethasone` sample

This is a **paired / blocked design**: cell line is a nuisance factor, not a replicate group. A correct model includes it (e.g. `~ cellLine + treatment` in DESeq2, or cell line as a blocking factor in limma/edgeR). Treating the 8 samples as 4-vs-4 unpaired throws away the strongest source of variance and is the most likely analysis mistake here.

`data/rsem.merged.gene_counts.tsv` — 58,735 genes × 8 samples of RSEM expected counts:

- First two columns are `gene_id` and `transcript_id(s)`; the remaining 8 are sample columns. Both must be dropped or set as identifiers before building a numeric matrix, and the `transcript_id(s)` header contains parentheses that `read.table`/`read_csv` will mangle unless name repair is disabled.
- Sample column order matches the row order of `experiment_table.csv`, but **verify by name rather than relying on position** when constructing the coldata/metadata object — a silent mismatch here produces plausible-looking but wrong results.
- Counts are **fractional** (RSEM expected counts, e.g. `82.33`). DESeq2 requires integers, so round before `DESeqDataSetFromMatrix`; the intended path for RSEM input is otherwise `tximport`, which is not applicable to this pre-merged gene-level file.
- Gene IDs are **unversioned Ensembl** (`ENSG00000000003`, no `.N` suffix) — no version stripping needed when mapping to symbols.
- 28,344 of 58,735 genes (~48%) are zero across all samples; the file is the full annotation, so pre-filtering low-count genes is expected.
