#!/usr/bin/env Rscript
#
# One-shot dependency install for airway_deseq2.qmd.
#
#   Rscript setup/install_packages.R
#
# Safe to re-run: already-installed packages are skipped. Everything installed
# here works offline afterwards -- the analysis document makes no network calls.

required <- c(
  # Core analysis
  "DESeq2",        # differential expression
  "apeglm",        # log2 fold change shrinkage (falls back to "normal" if absent)

  # Annotation -- optional at render time, the document degrades to Ensembl IDs
  "AnnotationDbi",
  "org.Hs.eg.db",  # human Ensembl ID -> gene symbol, local database

  # Plotting
  "ggplot2",
  "ggrepel",       # non-overlapping volcano labels (optional)
  "pheatmap",
  "RColorBrewer",

  # Pulled in by DESeq2, listed for clarity
  "matrixStats"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager from CRAN ...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All ", length(required), " packages already installed. Nothing to do.")
} else {
  message("Installing ", length(missing), " package(s): ",
          paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

# Report the final state so a partial install is obvious rather than silent.
status <- vapply(required, requireNamespace, logical(1), quietly = TRUE)

cat("\n--- dependency status ---\n")
for (pkg in required) {
  cat(sprintf("  %-16s %s\n", pkg,
              if (status[[pkg]]) as.character(packageVersion(pkg)) else "MISSING"))
}

cat("\nR version: ", R.version.string, "\n", sep = "")
if (requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Bioconductor: ", as.character(BiocManager::version()), "\n", sep = "")
}

if (any(!status)) {
  cat("\nStill missing: ", paste(required[!status], collapse = ", "), "\n", sep = "")
  cat("DESeq2 is required. The others are optional -- the document renders without\n")
  cat("them, using Ensembl IDs instead of symbols and 'normal' LFC shrinkage.\n")
  if (!status[["DESeq2"]]) quit(status = 1)
} else {
  cat("\nReady. Render with:  quarto render airway_deseq2.qmd\n")
}
