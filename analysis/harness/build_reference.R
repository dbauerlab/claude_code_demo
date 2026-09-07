#!/usr/bin/env Rscript
# Gold-standard reference: canonical DESeq2 workflow on the *published* Bioconductor
# `airway` dataset (Himes et al. 2014). This is deliberately independent of this repo's
# RSEM matrix so that "correctness" is judged against the community-standard result, not
# against any one branch. Because airway's counts come from a different quantification
# (Ensembl/HTSeq in the airway package vs RSEM here), we treat RANK / OVERLAP / DIRECTION /
# marker recovery as primary and absolute LFC magnitude as secondary.
suppressPackageStartupMessages({
  library(DESeq2); library(airway); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(1)
outdir <- file.path(dirname(dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1])))), "data")
if (!dir.exists(outdir)) outdir <- "analysis/data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data("airway")
se <- airway
# dex: 'untrt' (Untreated) vs 'trt' (Dexamethasone). Set Untreated as reference so the
# reported effect is Dexamethasone-vs-Untreated, matching every branch's contrast.
se$dex <- relevel(se$dex, ref = "untrt")
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- dds[rowSums(counts(dds)) >= 10, ]            # canonical vignette pre-filter
dds <- DESeq(dds)
res <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
res_sh <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm", res = res)

df <- as.data.frame(res)
df$gene_id <- sub("\\.\\d+$", "", rownames(df))     # strip any version suffix -> match repo IDs
df$log2FoldChange_shrunk <- res_sh$log2FoldChange[match(rownames(df), rownames(res_sh))]
df$symbol <- mapIds(org.Hs.eg.db, keys = df$gene_id, column = "SYMBOL",
                    keytype = "ENSEMBL", multiVals = "first")
df <- df[, c("gene_id","symbol","baseMean","log2FoldChange","log2FoldChange_shrunk",
             "lfcSE","stat","pvalue","padj")]
df <- df[order(df$padj), ]
write.csv(df, file.path(outdir, "reference_de.csv"), row.names = FALSE)

# Marker sanity check (Himes et al. known glucocorticoid responders)
markers <- c(DUSP1="ENSG00000120129", KLF15="ENSG00000163884", PER1="ENSG00000179094",
             CRISPLD2="ENSG00000103196", FKBP5="ENSG00000096060", SPARCL1="ENSG00000152583",
             CCL2="ENSG00000108691", ZBTB16="ENSG00000109906")
mk <- df[match(markers, df$gene_id), c("gene_id","symbol","log2FoldChange","padj")]
mk$expected <- names(markers)
cat("\n=== reference marker recovery (expect DUSP1/KLF15/SPARCL1/PER1 up, CCL2 down) ===\n")
print(mk, row.names = FALSE)
cat(sprintf("\nreference: %d genes tested, %d sig at padj<0.05\n",
            sum(!is.na(df$padj)), sum(df$padj < 0.05, na.rm = TRUE)))
writeLines(capture.output(sessionInfo()), file.path(outdir, "reference_sessionInfo.txt"))
cat("WROTE", file.path(outdir, "reference_de.csv"), "\n")
