#!/usr/bin/env Rscript
# Layer 3 quantitative comparison. Reads the harmonized per-branch DE tables + the airway
# reference, and computes: per-branch concordance vs reference (LFC correlation, top-N
# overlap, direction agreement, marker recovery), pairwise LFC correlation + significant-set
# Jaccard across branches, and the cross-branch consensus gene set. Writes tidy CSVs to
# analysis/data/ that the report and figures consume (no numbers are hand-entered downstream).
suppressWarnings(suppressMessages({library(stats)}))
D <- "analysis/data"; H <- file.path(D, "de_harmonized")
ref <- read.csv(file.path(D, "reference_de.csv"), stringsAsFactors = FALSE)
ref$gene_id <- sub("\\.\\d+$", "", ref$gene_id)
branches <- sub("\\.csv$", "", list.files(H, pattern = "\\.csv$"))
tabs <- setNames(lapply(branches, function(b) {
  x <- read.csv(file.path(H, paste0(b, ".csv")), stringsAsFactors = FALSE)
  x$gene_id <- sub("\\.\\d+$", "", x$gene_id); x[!duplicated(x$gene_id), ]
}), branches)

sig <- function(x) x$gene_id[which(x$padj < 0.05)]
jacc <- function(a, b) { u <- length(union(a,b)); if(!u) return(NA); length(intersect(a,b))/u }
ref_sig <- sig(ref)
markers <- c(SPARCL1="ENSG00000152583", DUSP1="ENSG00000120129", KLF15="ENSG00000163884",
             PER1="ENSG00000179094", CRISPLD2="ENSG00000103196", FKBP5="ENSG00000096060",
             ZBTB16="ENSG00000109906", CACNB2="ENSG00000165995", SAMHD1="ENSG00000101347",
             CCL2="ENSG00000108691")
up_markers <- setdiff(names(markers), "CCL2")  # CCL2 is the repressed control

## ---- per-branch vs reference ----
per <- do.call(rbind, lapply(branches, function(b) {
  x <- tabs[[b]]
  m <- merge(x[, c("gene_id","log2FoldChange","padj")],
             ref[, c("gene_id","log2FoldChange","padj")], by = "gene_id", suffixes = c("",".ref"))
  ok <- is.finite(m$log2FoldChange) & is.finite(m$log2FoldChange.ref)
  bs <- sig(x)
  # top-100 by padj overlap
  t100b <- head(x$gene_id[order(x$padj)], 100); t100r <- head(ref$gene_id[order(ref$padj)], 100)
  both_sig <- intersect(bs, ref_sig)
  ms <- merge(x[x$gene_id %in% both_sig, c("gene_id","log2FoldChange")],
              ref[, c("gene_id","log2FoldChange")], by = "gene_id", suffixes = c("",".ref"))
  data.frame(
    branch = b,
    n_tested = sum(!is.na(x$padj)),
    n_sig = length(bs),
    n_sig_lfc1 = sum(x$padj < 0.05 & abs(x$log2FoldChange) > 1, na.rm = TRUE),
    pearson_lfc = cor(m$log2FoldChange[ok], m$log2FoldChange.ref[ok]),
    spearman_lfc = cor(m$log2FoldChange[ok], m$log2FoldChange.ref[ok], method = "spearman"),
    top100_overlap = length(intersect(t100b, t100r)) / 100,
    jaccard_sig_ref = jacc(bs, ref_sig),
    dir_agree = mean(sign(ms$log2FoldChange) == sign(ms$log2FoldChange.ref)),
    markers_up_ok = mean(sapply(up_markers, function(g){
      r <- x[x$gene_id==markers[g],]; nrow(r)==1 && r$padj<0.05 && r$log2FoldChange>0 })),
    stringsAsFactors = FALSE)
}))
write.csv(per, file.path(D, "output_metrics.csv"), row.names = FALSE)

## ---- pairwise LFC (spearman) + sig Jaccard ----
all_b <- branches
S <- matrix(NA, length(all_b), length(all_b), dimnames = list(all_b, all_b))
J <- S
for (i in all_b) for (j in all_b) {
  m <- merge(tabs[[i]][,c("gene_id","log2FoldChange")], tabs[[j]][,c("gene_id","log2FoldChange")],
             by="gene_id"); ok <- is.finite(m[,2]) & is.finite(m[,3])
  S[i,j] <- cor(m[ok,2], m[ok,3], method="spearman")
  J[i,j] <- jacc(sig(tabs[[i]]), sig(tabs[[j]]))
}
write.csv(round(S,4), file.path(D,"pairwise_lfc_spearman.csv"))
write.csv(round(J,4), file.path(D,"pairwise_sig_jaccard.csv"))

## ---- consensus across branches ----
sig_lists <- lapply(tabs, sig)
tallies <- table(unlist(sig_lists))
consensus <- data.frame(n_branches = as.integer(tallies))
core <- names(tallies)[tallies == length(all_b)]                 # sig in ALL branches
core_and_ref <- intersect(core, ref_sig)
writeLines(c(sprintf("branches_compared_quantitatively: %d (%s)", length(all_b), paste(all_b, collapse=", ")),
             sprintf("reference_sig: %d", length(ref_sig)),
             sprintf("core_sig_in_all_branches: %d", length(core)),
             sprintf("core_also_sig_in_reference: %d (%.1f%%)", length(core_and_ref), 100*length(core_and_ref)/length(core)),
             sprintf("median_pairwise_spearman_LFC: %.3f", median(S[upper.tri(S)])),
             sprintf("median_pairwise_sig_jaccard: %.3f", median(J[upper.tri(J)]))),
           file.path(D, "consensus_summary.txt"))

## ---- marker recovery table (all branches + reference) ----
mrec <- do.call(rbind, lapply(c(list(reference=ref), tabs), function(x) NULL))
mk <- do.call(rbind, lapply(names(markers), function(g){
  row <- data.frame(marker=g, gene_id=markers[g], stringsAsFactors=FALSE)
  for (src in c("reference", all_b)) {
    tb <- if (src=="reference") ref else tabs[[src]]
    r <- tb[tb$gene_id==markers[g],]
    row[[paste0("lfc_",src)]] <- if(nrow(r)==1) round(r$log2FoldChange,2) else NA
    row[[paste0("padj_",src)]] <- if(nrow(r)==1) r$padj else NA
  }
  row
}))
write.csv(mk, file.path(D,"marker_recovery.csv"), row.names=FALSE)

cat("=== per-branch vs airway reference ===\n"); print(per, digits=3, row.names=FALSE)
cat("\n"); writeLines(readLines(file.path(D,"consensus_summary.txt")))
cat("\nWROTE output_metrics.csv, pairwise_*.csv, marker_recovery.csv, consensus_summary.txt\n")
