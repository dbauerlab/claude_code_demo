# claude_code_demo — airway RNA-seq differential expression

Bulk RNA-seq differential expression demo on the **airway** dataset (Himes et al. — human airway
smooth muscle cells, dexamethasone vs. untreated), quantified with RSEM and merged by
nf-core/rnaseq.

Two implementations of the same analysis are provided — an R/Quarto document and a Python
script. They share the dataset, the paired design, and every gotchas-handling step; pick
whichever you prefer to run.

## Contents

| Path | Description |
|---|---|
| `airway_deseq2.qmd` | **R path.** Quarto document; runs the whole workflow and writes all outputs. |
| `analyze.py` | **Python path.** Standalone script (PyDESeq2); same workflow, run with `python analyze.py`. |
| `data/experiment_table.csv` | Sample sheet — 8 samples, `ID, cellLine, treatment`. The only place the design lives. |
| `data/rsem.merged.gene_counts.tsv` | 58,735 genes × 8 samples of RSEM expected counts. |
| `setup/install_packages.R` | One-shot R dependency install. |
| `setup/requirements.txt` | Python dependencies. |

Either path creates `results/` (TSV tables) and `figures/` (PNG + PDF per figure). Both are
gitignored. For the R path the rendered `airway_deseq2.html` is the shareable product.

## Run it — R

Requires R (≥ 4.3 recommended) and [Quarto](https://quarto.org/docs/get-started/).

```bash
Rscript setup/install_packages.R
```
```bash
quarto render airway_deseq2.qmd
```

DESeq2 is the only hard requirement. `org.Hs.eg.db` (gene symbols), `apeglm` (fold-change
shrinkage) and `ggrepel` (plot labels) are optional — without them the document still renders,
falling back to Ensembl IDs, `normal` shrinkage, and plain labels respectively.

## Run it — Python

Requires Python (≥ 3.9). [PyDESeq2](https://pydeseq2.readthedocs.io/) is a faithful Python
reimplementation of DESeq2; results track the R analysis closely.

```bash
pip install -r setup/requirements.txt
```
```bash
python analyze.py
```

PyDESeq2 is the only hard requirement. Gene symbols are optional and, unlike R, come from
`pyensembl` (a one-time offline data install) rather than a bundled database — without it the
script labels genes with Ensembl IDs:

```bash
pip install pyensembl && pyensembl install --release 110 --species homo_sapiens
```

> **Note:** CLAUDE.md specifies R/DESeq2 as the project's analysis path. The Python script is an
> alternative kept on this `demo/*` branch; the R Quarto document remains the
> convention-compliant deliverable.

## The analysis

The design is **paired**: 4 cell lines × (`Untreated`, `Dexamethasone`), one sample each. Cell
line is a nuisance factor, not a replicate group, so the model blocks on it:

```r
~ cellLine + treatment      # treatment last, so it is the coefficient of interest
```

`Untreated` is set explicitly as the reference level, which means **positive log2 fold change =
induced by dexamethasone**.

The document validates its own inputs with `stopifnot()` assertions rather than printed notes —
sample/metadata alignment is checked by name, the reference level is asserted, and the gene count
after filtering is bounded. Any of these failing stops the render instead of producing
plausible-looking but wrong results. A final diagnostic section refits `~ treatment` alone to
quantify what discarding the pairing would cost.

## Workshop note

Each participant works on their own `demo/<name>` branch off `main`. Keep new work on your
`demo/*` branch.
