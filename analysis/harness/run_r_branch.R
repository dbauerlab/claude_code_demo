#!/usr/bin/env Rscript
# Faithfully re-run one R/Quarto branch and harvest its DE results table WITHOUT altering
# its method. We tangle the .qmd to R (knitr::purl), execute it in the branch's own dir
# (each carries its own data/), then introspect the environment for the results object the
# branch itself produced. Nothing about the statistics is reimplemented here — we only read
# back whatever the branch computed. Errors in late plotting/annotation chunks are tolerated:
# objects created before the failure still persist and are harvested.
args <- commandArgs(TRUE)
branch_dir <- normalizePath(args[1]); qmd <- args[2]; out_csv <- normalizePath(args[3], mustWork = FALSE)
log <- function(...) cat(sprintf(...), "\n")
suppressWarnings(suppressMessages(library(knitr)))

owd <- setwd(branch_dir)
tmpR <- tempfile(fileext = ".R")
purl(qmd, output = tmpR, quiet = TRUE, documentation = 0L)
# headless graphics: swallow all plotting so device calls never abort the run
pdf(file = nullfile()); on.exit(try(dev.off(), silent = TRUE), add = TRUE)
e <- new.env()
err <- tryCatch({ sys.source(tmpR, envir = e); NA_character_ },
                error = function(x) conditionMessage(x))
setwd(owd)
if (!is.na(err)) log("  (script stopped in a later chunk: %s)", substr(err, 1, 80))

# Find every data-frame-like object carrying a log2FoldChange + padj; pick the fullest table.
cand <- list()
for (nm in ls(e)) {
  obj <- tryCatch(get(nm, e), error = function(x) NULL)
  if (!is.data.frame(obj) && !inherits(obj, c("DataFrame","DFrame","data.table","tbl_df"))) next
  df <- tryCatch(as.data.frame(obj), error = function(x) NULL)
  if (is.null(df) || is.null(nrow(df)) || nrow(df) == 0L) next
  cn <- tolower(names(df))
  if (any(grepl("log2foldchange|log2fc", cn)) && any(grepl("padj|adj.*p|qvalue", cn)))
    cand[[nm]] <- df
}
if (!length(cand)) stop("no DE results object found in ", branch_dir)
# Prefer the object with the most rows (the full per-gene table, not a top-N slice).
best <- names(cand)[which.max(vapply(cand, nrow, 1L))]
df <- cand[[best]]
log("  harvested object '%s' (%d rows) from {%s}", best, nrow(df), paste(names(cand), collapse=", "))

pick <- function(df, pats) { cn <- tolower(names(df)); for (p in pats) { h <- which(grepl(p, cn)); if (length(h)) return(names(df)[h[1]]) }; NA }
gid <- pick(df, c("^gene_id$","^gene$","^ensembl","^id$"))
gene_id <- if (!is.na(gid)) as.character(df[[gid]]) else rownames(df)
lfc  <- pick(df, c("log2foldchange_shrunk"))      # prefer shrunken if present, else raw
lfc_raw <- pick(df, c("^log2foldchange$","log2foldchange(?!_)","log2fc"))
out <- data.frame(
  gene_id = sub("\\.\\d+$", "", gene_id),
  baseMean = if (!is.na(pick(df,c("basemean")))) df[[pick(df,c("basemean"))]] else NA_real_,
  log2FoldChange = df[[if (!is.na(lfc_raw)) lfc_raw else lfc]],
  log2FoldChange_shrunk = if (!is.na(lfc)) df[[lfc]] else NA_real_,
  pvalue = if (!is.na(pick(df,c("^pvalue$","p_value","pval")))) df[[pick(df,c("^pvalue$","p_value","pval"))]] else NA_real_,
  padj = df[[pick(df, c("^padj$","adj.*p","qvalue"))]],
  stringsAsFactors = FALSE)
out <- out[!is.na(out$gene_id) & out$gene_id != "", ]
out <- out[!duplicated(out$gene_id), ]
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(out, out_csv, row.names = FALSE)
log("  WROTE %s (%d genes, %d sig padj<0.05)", out_csv, nrow(out), sum(out$padj < 0.05, na.rm = TRUE))
