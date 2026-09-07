#!/usr/bin/env Rscript
# Layer 1: quantitative + rubric scoring of each branch's CLAUDE.md context file.
# Rubric = presence (case-insensitive regex) of key correct/important analysis facts an
# agent should have captured. Objective keyword detection; the qualitative read lives in the
# report prose. Also emits size/structure descriptors. Scored on the materialized branches.
W <- "analysis/work/branches"; D <- "analysis/data"
branches <- c("Catarina","Olivia","santosnuno","okp-23","becky","tiph","barretj","pateld1","davidlvb","marga")
rubric <- c(
  dataset_airway      = "airway|himes",
  paired_design       = "paired|blocking|block|nuisance",
  formula             = "celll?ine *\\+ *treatment|cell *\\+ *dex",
  ref_untreated       = "reference.{0,20}untreated|untreated.{0,20}reference|ref *= *.?untreated|relevel",
  signflip_warning    = "alphabetical|sign of|flip|invert|silently",
  round_counts        = "round",
  round_rationale     = "multi-?mapping|floor|trunc|as\\.integer|probabilis|expected count",
  prefilter           = "pre-?filter|filter|rowsums|low-?count|all-?zero",
  join_by_id          = "join|match.{0,10}id|by name|not.{0,6}position|verify by",
  drop_transcript     = "transcript_id|name.?repair|drop.{0,15}transcript|second column",
  ensembl_unversioned = "unversioned|without.{0,6}version|no.{0,6}version|\\.n suffix|no `?\\.",
  symbol_mapping      = "org\\.hs\\.eg\\.db|biomart|gene symbol|annotables|map.{0,10}symbol",
  libsize_var         = "library size|lib.{0,4}size|size factor|normali",
  marker_genes        = "dusp1|klf15|crispld2|per1|sparcl1|ccl2|fkbp5",
  lang_convention     = "quarto|\\.qmd|deseq2|pydeseq2|deliverable|language",
  quant_facts         = "58,?735|28,?344|11,?117|14\\.9|30\\.3|20,?783|20,?803"
)
rows <- lapply(branches, function(b) {
  f <- file.path(W, b, "CLAUDE.md"); txt <- tolower(paste(readLines(f, warn = FALSE), collapse = "\n"))
  raw <- paste(readLines(f, warn = FALSE), collapse = "\n")
  sc <- sapply(rubric, function(p) as.integer(grepl(p, txt, perl = TRUE)))
  data.frame(branch = b,
    words = lengths(gregexpr("\\S+", raw)),
    lines = length(readLines(f, warn = FALSE)),
    headings = length(gregexpr("(^|\n)#+ ", raw)[[1]]),
    code_blocks = length(gregexpr("```", raw)[[1]]) %/% 2,
    coverage = sum(sc), t(sc), stringsAsFactors = FALSE, check.names = FALSE)
})
tab <- do.call(rbind, rows)
write.csv(tab, file.path(D, "claudemd_rubric.csv"), row.names = FALSE)
cat("=== CLAUDE.md descriptors + rubric coverage (/16) ===\n")
print(tab[, c("branch","words","lines","headings","code_blocks","coverage")], row.names = FALSE)
cat("\nrubric items:", paste(names(rubric), collapse=", "), "\n")
cat("mean coverage:", round(mean(tab$coverage),1), "/ 16 ; range", min(tab$coverage), "-", max(tab$coverage), "\n")
