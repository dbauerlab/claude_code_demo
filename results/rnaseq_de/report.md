# ClawBio RNA-seq Differential Expression Report

**Date**: 2026-08-05 16:20 UTC
**Samples**: 8
**Genes (pre-filter)**: 58735
**Genes (post-filter)**: 16477
**Formula**: `~ cellLine + treatment`
**Contrast**: `treatment,Dexamethasone,Untreated`
**Backend used**: `pydeseq2`
**LFC shrinkage**: `applied`
**LFC shrinkage coefficient**: `treatment[T.Dexamethasone]`


## Pre-DE QC + PCA

- QC summary: `tables/qc_summary.csv`
- PCA figure: `figures/pca.png`

## Differential Expression

- Full results: `tables/de_results.csv`
- Volcano plot: `figures/volcano.png`
- MA plot: `figures/ma_plot.png`

### Top Genes (by adjusted p-value)

| Gene | log2FoldChange | padj |
|---|---:|---:|
| ENSG00000152583 | 4.563 | 5.490e-168 |
| ENSG00000101347 | 3.746 | 8.697e-165 |
| ENSG00000211445 | 3.714 | 3.403e-147 |
| ENSG00000189221 | 3.333 | 1.526e-145 |
| ENSG00000120129 | 2.938 | 4.011e-143 |
| ENSG00000196136 | 3.202 | 9.781e-131 |
| ENSG00000165995 | 3.304 | 1.094e-115 |
| ENSG00000154734 | 2.331 | 1.812e-101 |
| ENSG00000162614 | 2.016 | 3.051e-100 |
| ENSG00000163884 | 4.432 | 7.731e-96 |

## Reproducibility

- Commands: `reproducibility/commands.sh`
- Environment: `reproducibility/environment.yml`
- Checksums: `reproducibility/checksums.sha256`

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
