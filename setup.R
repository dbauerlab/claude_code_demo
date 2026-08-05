# setup.R — one-time package installation for the dexamethasone DE analysis.
#
# Idempotent: each package is checked before installing, so re-running this is cheap.
# Run once with:  Rscript setup.R
#
# Expect ~10-20 minutes on a cold install. DESeq2 pulls a large dependency tree and
# org.Hs.eg.db is a sizeable annotation database.

# Set a CRAN mirror explicitly. Without this, a non-interactive `Rscript setup.R` can fail with
# "trying to use CRAN without setting a mirror" when resolving CRAN dependencies.
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

needed <- c(
  # Bioconductor
  "DESeq2",         # differential expression
  "apeglm",         # backs lfcShrink(type = "apeglm")
  "org.Hs.eg.db",   # offline human Ensembl -> symbol annotation
  "AnnotationDbi",  # mapIds() interface to the above
  # CRAN
  "knitr",          # Quarto's R engine -- required to execute the .qmd chunks
  "rmarkdown",      # ditto; Quarto shells out to it for R rendering
  "ggplot2",
  "ggrepel",        # non-overlapping gene labels on the volcano plot
  "pheatmap",
  "RColorBrewer"
)

missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
} else {
  message("All packages already installed.")
}

# Fail loudly if anything did not install, rather than letting the .qmd fail mid-render.
still_missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}

message("Setup complete. Render with: quarto render dexamethasone_de.qmd")
