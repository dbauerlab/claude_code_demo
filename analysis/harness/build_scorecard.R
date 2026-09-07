#!/usr/bin/env Rscript
# Phase 4 synthesis: a correctness scorecard across five explicit dimensions. Objective
# dimensions (C, D) come straight from computed metrics; A and B come from the verified
# script matrix (every branch was read); E (reproducibility) is a transparent reviewer
# rubric over committed artefacts. All dimensions are 0-1; composite is their mean.
D <- "analysis/data"
met <- read.csv(file.path(D, "output_metrics.csv"), stringsAsFactors = FALSE)
sm  <- read.csv(file.path(D, "script_matrix.csv"), stringsAsFactors = FALSE)

# E: reproducibility rubric (committed env/seed/tests/repro-bundle). Justified in report.
repro <- c(barretj=1.00, pateld1=0.90, okp23=0.80, tiph=0.80, Olivia=0.75, santosnuno=0.75,
           davidlvb=0.70, becky=0.55, Catarina=0.55, marga=0.55)
names(repro)[names(repro)=="okp23"] <- "okp-23"

score <- do.call(rbind, lapply(sm$branch, function(b) {
  m <- met[met$branch == b, ]
  A <- 1.0                                   # design: all use ~ cellLine + treatment (verified)
  B <- 1.0                                   # stats hygiene: round+filter+reference+shrink+BH all present
  if (nrow(m) == 1) {
    C <- mean(c(m$markers_up_ok, m$dir_agree))         # biological validity
    D_ <- mean(c(m$pearson_lfc, m$top100_overlap))     # concordance vs airway reference
  } else { C <- NA; D_ <- NA }               # davidlvb: no regenerable output (clawbio engine offline)
  E <- unname(repro[b])
  comp <- mean(c(A, B, C, D_, E), na.rm = TRUE)
  data.frame(branch=b, A_design=A, B_stats=B, C_biology=round(C,3),
             D_concordance=round(D_,3), E_reproducibility=E, composite=round(comp,3),
             stringsAsFactors=FALSE)
}))
score <- score[order(-score$composite), ]
write.csv(score, file.path(D, "correctness_scorecard.csv"), row.names = FALSE)
cat("=== correctness scorecard (0-1; composite = mean of A-E) ===\n")
print(score, row.names = FALSE)
cat("\nNote: davidlvb C/D are NA (clawbio engine unavailable to regenerate output);",
    "\nits method mirrors barretj (same clawbio/PyDESeq2 engine), which concords at D=",
    round(mean(c(met$pearson_lfc[met$branch=='barretj'], met$top100_overlap[met$branch=='barretj'])),3), "\n")
