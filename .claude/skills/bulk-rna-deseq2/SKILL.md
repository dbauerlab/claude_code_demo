---
name: bulk-rna-deseq2
description: >
  Use when performing bulk RNA-seq differential expression in R with DESeq2 from a
  gene-level counts matrix and a sample sheet. Covers reading counts (including
  non-integer RSEM/tximport counts), building a DESeqDataSet, running the standard
  DESeq2 workflow, and producing QC, results tables, and figures. Deliverables are
  reproducible Quarto (.qmd) documents.
---

# Bulk RNA-seq Differential Expression (R / DESeq2)

## When to use
A gene-level counts matrix + a sample sheet, and the user wants differential expression.

## Standard workflow
1. **Load inputs.** Read the counts matrix; drop non-count columns (e.g. `transcript_id(s)`);
   set gene IDs as rownames. Read the sample sheet.
2. **Handle non-integer counts.** RSEM/tximport counts are estimates — `round()` the matrix
   to integers before building the DESeqDataSet (or import via `tximport`).
3. **Align samples.** Reorder counts columns to match the sample-sheet row order; assert they
   are identical before proceeding.
4. **Factors & reference levels.** Convert design variables to factors. Explicitly set the
   control/reference level with `relevel()` (e.g. `Untreated`).
5. **Build & run.** `DESeqDataSetFromMatrix(countData, colData, design)` → `DESeq()`.
   Respect a paired/blocked design if donors/batches exist (e.g. `~ subject + condition`).
6. **QC.** Pre-filter low-count genes. `vst()`/`rlog` transform, then PCA and a
   sample-distance heatmap to check clustering by condition.
7. **Results.** `results()` with a stated alpha (e.g. 0.05); `lfcShrink()` for ranking/plots.
   Order by adjusted p-value; report how many genes pass FDR.
8. **Figures.** MA plot, volcano plot, and a heatmap of the top variable/DE genes.
9. **Export.** Save the full results table as CSV and figures under `results/`.
10. **Reproducibility.** End with `sessionInfo()`.

## Output contract
Produce a single Quarto (`.qmd`) document that renders end-to-end with no manual steps,
with clear section headings and a short interpretation under each figure.

## Common pitfalls
- Forgetting to round non-integer counts → DESeq2 errors.
- Reference level defaults to alphabetical → wrong sign on log2FC. Always `relevel()`.
- Counts columns not matched to the sample sheet → silently wrong results.
