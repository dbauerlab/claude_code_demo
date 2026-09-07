#!/usr/bin/env Rscript
# Generate every figure the report embeds, from the computed CSVs + harmonized tables.
suppressWarnings(suppressMessages({library(ggplot2); library(reshape2)}))
D <- "analysis/data"; H <- file.path(D,"de_harmonized"); F <- "analysis/figures"
dir.create(F, showWarnings = FALSE, recursive = TRUE)
ok <- "#2a7ab0"; hot <- "#c1442e"; mid <- "#f0f0f0"
th <- theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))
sv <- function(p, f, w=8, h=5.2) ggsave(file.path(F,f), p, width=w, height=h, dpi=130, bg="white")

## 1. CLAUDE.md rubric coverage heatmap
r <- read.csv(file.path(D,"claudemd_rubric.csv"), check.names=FALSE)
items <- r[, 7:ncol(r)]
m <- melt(cbind(branch=r$branch, items), id.vars="branch")
m$branch <- factor(m$branch, levels=r$branch[order(-r$coverage)])
p1 <- ggplot(m, aes(variable, branch, fill=factor(value))) +
  geom_tile(color="white", linewidth=0.6) +
  scale_fill_manual(values=c("0"=mid,"1"=ok), labels=c("absent","present"), name=NULL) +
  labs(title="Layer 1 — CLAUDE.md context coverage",
       subtitle="Presence of each key analysis fact (rows ordered by total coverage /16)",
       x=NULL, y=NULL) + th
sv(p1, "fig1_claudemd_coverage.png", 10, 5.4)

## 2. CLAUDE.md size vs coverage
r$branch <- factor(r$branch, levels=r$branch[order(r$coverage)])
p2 <- ggplot(r, aes(coverage, branch)) +
  geom_col(aes(fill=words), width=0.7) +
  geom_text(aes(label=paste0(coverage,"/16")), hjust=-0.15, size=3.4) +
  scale_fill_gradient(low="#cfe3f2", high=ok, name="words") +
  scale_x_continuous(limits=c(0,17), expand=expansion(mult=c(0,0.02))) +
  labs(title="Layer 1 — context richness", x="rubric coverage (/16)", y=NULL) + th +
  theme(axis.text.x=element_text(angle=0))
sv(p2, "fig2_claudemd_size.png", 7.5, 4.6)

## 3. pairwise LFC Spearman heatmap
S <- as.matrix(read.csv(file.path(D,"pairwise_lfc_spearman.csv"), row.names=1, check.names=FALSE))
mS <- melt(S);
p3 <- ggplot(mS, aes(Var1, Var2, fill=value)) + geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f",value)), size=2.7) +
  scale_fill_gradient2(low="white", mid="#cfe3f2", high=ok, midpoint=0.9, limits=c(0.8,1), name="Spearman") +
  labs(title="Layer 3 — pairwise log2FC rank concordance", subtitle="Spearman correlation of per-gene LFC between branches", x=NULL, y=NULL) + th
sv(p3, "fig3_pairwise_lfc.png", 8.2, 6.4)

## 4. significant-set Jaccard heatmap
J <- as.matrix(read.csv(file.path(D,"pairwise_sig_jaccard.csv"), row.names=1, check.names=FALSE))
mJ <- melt(J)
p4 <- ggplot(mJ, aes(Var1, Var2, fill=value)) + geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f",value)), size=2.7) +
  scale_fill_gradient(low="white", high="#2e8b57", limits=c(0.7,1), name="Jaccard") +
  labs(title="Layer 3 — significant-gene set agreement", subtitle="Jaccard overlap of padj<0.05 gene sets (DESeq2 vs PyDESeq2 vs clawbio)", x=NULL, y=NULL) + th
sv(p4, "fig4_sig_jaccard.png", 8.2, 6.4)

## 5. concordance vs airway reference
met <- read.csv(file.path(D,"output_metrics.csv"))
mm <- melt(met[,c("branch","pearson_lfc","spearman_lfc","top100_overlap","jaccard_sig_ref")], id.vars="branch")
lab <- c(pearson_lfc="LFC Pearson", spearman_lfc="LFC Spearman", top100_overlap="top-100 overlap", jaccard_sig_ref="sig Jaccard")
mm$variable <- lab[as.character(mm$variable)]
mm$branch <- factor(mm$branch, levels=met$branch[order(met$pearson_lfc)])
p5 <- ggplot(mm, aes(value, branch, fill=variable)) +
  geom_col(position=position_dodge(0.75), width=0.7) +
  scale_fill_manual(values=c("LFC Pearson"=ok,"LFC Spearman"="#7fb2d6","top-100 overlap"="#e08a3c","sig Jaccard"="#2e8b57"), name=NULL) +
  scale_x_continuous(limits=c(0,1), expand=expansion(mult=c(0,0.02))) +
  labs(title="Layer 3 — concordance with published airway reference", x="score (0-1)", y=NULL) + th +
  theme(axis.text.x=element_text(angle=0))
sv(p5, "fig5_vs_reference.png", 8.5, 5.2)

## 6. N genes tested vs significant
ng <- melt(met[,c("branch","n_tested","n_sig")], id.vars="branch")
ng$branch <- factor(ng$branch, levels=met$branch[order(met$n_tested)])
p6 <- ggplot(ng, aes(value, branch, fill=variable)) +
  geom_col(position=position_dodge(0.7), width=0.68) +
  scale_fill_manual(values=c(n_tested="#b9c6cf", n_sig=hot), labels=c("genes tested","significant (padj<0.05)"), name=NULL) +
  labs(title="Layer 3 — genes tested vs significant", subtitle="Filter choice moves genes-tested a lot; significant count barely at all", x="genes", y=NULL) + th +
  theme(axis.text.x=element_text(angle=0))
sv(p6, "fig6_ngenes.png", 8.2, 5)

## 7. marker-gene recovery heatmap (LFC), branches + reference
mk <- read.csv(file.path(D,"marker_recovery.csv"), check.names=FALSE)
lfc_cols <- grep("^lfc_", names(mk), value=TRUE)
ml <- melt(mk[,c("marker",lfc_cols)], id.vars="marker", variable.name="src", value.name="lfc")
ml$src <- sub("^lfc_","",ml$src)
ml$src <- factor(ml$src, levels=c("reference", setdiff(unique(ml$src),"reference")))
ml$marker <- factor(ml$marker, levels=rev(mk$marker))
p7 <- ggplot(ml, aes(src, marker, fill=lfc)) + geom_tile(color="white") +
  geom_text(aes(label=ifelse(is.na(lfc),"NA",sprintf("%.1f",lfc))), size=2.9) +
  scale_fill_gradient2(low=ok, mid="white", high=hot, midpoint=0, na.value="grey85", name="log2FC") +
  labs(title="Biological validity — known glucocorticoid responders",
       subtitle="Dexamethasone should induce all but CCL2 (repressed). 'reference' = airway package.", x=NULL, y=NULL) + th
sv(p7, "fig7_markers.png", 9, 5.2)

## 8. LFC vs reference scatter (two representative branches: tightest + an engine port)
ref <- read.csv(file.path(D,"reference_de.csv")); ref$gene_id <- sub("\\.\\d+$","",ref$gene_id)
pick <- c("tiph","pateld1")
sc <- do.call(rbind, lapply(pick, function(b){
  x <- read.csv(file.path(H,paste0(b,".csv"))); x$gene_id <- sub("\\.\\d+$","",x$gene_id)
  z <- merge(x[,c("gene_id","log2FoldChange")], ref[,c("gene_id","log2FoldChange")], by="gene_id", suffixes=c("",".ref"))
  z$branch <- b; z[is.finite(z$log2FoldChange)&is.finite(z$log2FoldChange.ref),]
}))
p8 <- ggplot(sc, aes(log2FoldChange.ref, log2FoldChange)) +
  geom_abline(slope=1, intercept=0, color="grey60", linetype=2) +
  geom_point(alpha=0.12, size=0.5, color=ok) +
  facet_wrap(~branch) +
  labs(title="Layer 3 — per-gene log2FC vs airway reference",
       subtitle="tiph (R/DESeq2) and pateld1 (PyDESeq2); dashed line = identity",
       x="reference log2FC (airway package)", y="branch log2FC") +
  theme_minimal(base_size=12) + theme(strip.text=element_text(face="bold"))
sv(p8, "fig8_lfc_scatter.png", 8.5, 4.6)

## 9. correctness scorecard heatmap
sco <- read.csv(file.path(D,"correctness_scorecard.csv"))
sl <- melt(sco, id.vars="branch", variable.name="dim", value.name="v")
sl <- sl[sl$dim!="composite",]
sl$dim <- factor(sl$dim, levels=c("A_design","B_stats","C_biology","D_concordance","E_reproducibility"))
sl$branch <- factor(sl$branch, levels=rev(sco$branch))
p9 <- ggplot(sl, aes(dim, branch, fill=v)) + geom_tile(color="white") +
  geom_text(aes(label=ifelse(is.na(v),"n/a",sprintf("%.2f",v))), size=3) +
  scale_fill_gradient(low="#fbe6df", high="#2e8b57", limits=c(0.5,1), na.value="grey85", name="score") +
  labs(title="Correctness scorecard", subtitle="A design · B stats hygiene · C biology · D concordance vs ref · E reproducibility", x=NULL, y=NULL) + th
sv(p9, "fig9_scorecard.png", 8.2, 5)

cat("wrote", length(list.files(F, pattern="png$")), "figures to", F, "\n")
