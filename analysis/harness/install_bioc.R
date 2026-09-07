# Best-effort Bioconductor + CRAN install for the comparison harness + reference run.
options(warn = 1, Ncpus = max(1, parallel::detectCores() - 1))
cran <- c("BiocManager","pheatmap","UpSetR","ggrepel","matrixStats")
for (p in cran) if (!requireNamespace(p, quietly = TRUE))
  tryCatch(install.packages(p, repos = "https://cloud.r-project.org"),
           error = function(e) message("CRAN FAIL ", p, ": ", conditionMessage(e)))
bioc <- c("DESeq2","airway","org.Hs.eg.db","apeglm","ashr")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE))
  tryCatch(BiocManager::install(p, update = FALSE, ask = FALSE),
           error = function(e) message("BIOC FAIL ", p, ": ", conditionMessage(e)))
cat("\n=== INSTALL RESULTS ===\n")
for (p in c(cran, bioc))
  cat(sprintf("%-16s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("DONE\n")
