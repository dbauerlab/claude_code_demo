#!/usr/bin/env python
"""
Dexamethasone response in airway smooth muscle — bulk RNA-seq differential expression.

Python port of ``airway_deseq2.qmd``, using PyDESeq2 (a faithful reimplementation of
DESeq2) in place of the R/Bioconductor stack. Same dataset, same paired design, same
gotchas handled the same way.

Run it:

    python analyze.py

Writes TSV tables to ``results/`` and figures (PNG + PDF) to ``figures/``, and prints a
running set of self-checks. As in the R version, every consistency check is an assertion:
if an input changes shape, a sample is mismatched, or the reference level flips, the script
stops loudly instead of producing plausible-but-wrong numbers.

Design: ``~ cellLine + treatment`` — 4 cell lines each contribute one Untreated and one
Dexamethasone sample, so cell line is a nuisance factor blocked in the model, not a
replicate group. ``Untreated`` is the reference, so a positive log2 fold change means the
gene is induced by dexamethasone.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")  # headless: render straight to files, no display needed
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

RANDOM_SEED = 20260805

# Known dexamethasone responders in airway smooth muscle, used as a biological positive
# control. IDs were verified present in the counts file and passing the count filter
# before being hard-coded here — they are checked, not recalled — and each is still
# guarded at plot time in case a future counts file drops one.
MARKER_GENES = {
    "CRISPLD2": "ENSG00000103196",
    "DUSP1": "ENSG00000120129",
    "PER1": "ENSG00000179094",
    "KLF15": "ENSG00000163884",
    "FKBP5": "ENSG00000096060",
    "ZBTB16": "ENSG00000109906",
    "SPARCL1": "ENSG00000152583",
}

ID_COLS = ["gene_id", "transcript_id(s)"]
COEF_CONTRAST = ["treatment", "Dexamethasone", "Untreated"]
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------


def rule(title: str) -> None:
    """Print a section banner so the console log reads like the rendered document."""
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def savefig(fig, outdir: Path, name: str) -> None:
    """Write a figure as both PNG and PDF, mirroring the R document's dual output."""
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{name}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


def build_symbol_map(gene_ids, release: int):
    """
    Map unversioned Ensembl gene IDs to symbols via pyensembl, offline.

    Optional, exactly like org.Hs.eg.db in the R version: requires a one-time
    ``pyensembl install --release <N> --species homo_sapiens``. If pyensembl or its data
    are unavailable, returns an all-None map and the caller falls back to Ensembl IDs.
    """
    try:
        from pyensembl import EnsemblRelease

        data = EnsemblRelease(release=release, species="homo_sapiens")
        # Force a lookup so a missing local cache fails here, not later.
        _ = data.gene_by_id(gene_ids[0]) if len(gene_ids) else None

        mapping = {}
        for g in gene_ids:
            try:
                mapping[g] = data.gene_by_id(g).gene_name or None
            except Exception:
                mapping[g] = None
        return mapping, True
    except Exception as exc:  # pyensembl absent, or release data not installed
        print(f"  symbol mapping unavailable ({type(exc).__name__}) — using Ensembl IDs")
        print(f"  to enable: pip install pyensembl && pyensembl install "
              f"--release {release} --species homo_sapiens")
        return {g: None for g in gene_ids}, False


def make_labeller(symbol_map):
    """Return f(id)->symbol-or-id and f(id)->symbol-or-None helpers."""

    def label(gid):
        s = symbol_map.get(gid)
        return s if s else gid

    def symbol(gid):
        return symbol_map.get(gid)

    return label, symbol


def pca_2d(mat_samples_by_genes: np.ndarray, ntop: int = 500):
    """
    Replicate DESeq2::plotPCA: take the ``ntop`` most-variable genes, centre, and take the
    top-2 principal components via SVD. Returns (scores[:, :2], percent_variance[:2]).
    """
    variances = mat_samples_by_genes.var(axis=0)
    top_idx = np.argsort(variances)[::-1][:ntop]
    sub = mat_samples_by_genes[:, top_idx]
    centred = sub - sub.mean(axis=0)
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    scores = u * s
    var_explained = (s ** 2) / np.sum(s ** 2)
    return scores[:, :2], var_explained[:2]


# --------------------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------------------


def load_inputs(data_dir: Path):
    """Load the sample sheet and counts, and build a validated gene x sample matrix."""
    rule("1. Load and validate inputs")

    coldata = pd.read_csv(data_dir / "experiment_table.csv")
    # pandas preserves the 'transcript_id(s)' header verbatim — no name mangling — and
    # handles the file's Windows line endings transparently.
    raw = pd.read_csv(data_dir / "rsem.merged.gene_counts.tsv", sep="\t")

    print(f"sample sheet: {coldata.shape[0]} rows x {coldata.shape[1]} cols")
    print(f"counts file:  {raw.shape[0]} rows x {raw.shape[1]} cols")
    print(f"counts columns: {', '.join(raw.columns)}")

    assert all(c in raw.columns for c in ID_COLS), "missing identifier column(s)"
    assert raw["gene_id"].is_unique, "duplicate gene_id values"

    sample_ids = coldata["ID"].tolist()
    file_samples = [c for c in raw.columns if c not in ID_COLS]

    # Fail rather than guess if the sample sets differ at all.
    assert set(file_samples) == set(sample_ids), "counts columns != sample sheet IDs"

    # Match samples by NAME, not position: reorder the count columns to the sample sheet.
    # Order happens to agree here, but this is what keeps the script correct if the counts
    # file is ever regenerated with columns in a different order.
    counts_gxs = raw.set_index("gene_id")[sample_ids]
    assert list(counts_gxs.columns) == sample_ids
    print(f"sample columns aligned to sample sheet by name: True ({len(sample_ids)} samples)")

    return coldata, counts_gxs


def show_design(coldata: pd.DataFrame) -> None:
    """Print the design cross-tab; every cell must be 1 (the paired design)."""
    rule("Design")
    tab = pd.crosstab(coldata["cellLine"], coldata["treatment"])
    print(tab.to_string())
    assert (tab.values == 1).all(), "design is not a balanced one-per-cell paired design"
    print(f"\nbalanced paired design confirmed: {tab.shape[0]} cell lines x "
          f"{tab.shape[1]} treatments, one sample each")


def prepare_counts(counts_gxs: pd.DataFrame):
    """Round RSEM expected counts to integers and pre-filter low-count genes."""
    rule("2. Prepare counts")

    # RSEM reports fractional 'expected counts'. DESeq2/PyDESeq2 require integers. Round
    # with numpy round() — NOT int() truncation, which biases every fractional count down.
    n_nonint = int((counts_gxs != counts_gxs.round()).any(axis=1).sum())
    print(f"gene rows with a non-integer expected count: {n_nonint}")

    counts_int = counts_gxs.round().astype("int64")
    assert (counts_int.values >= 0).all(), "negative counts"

    n_before = counts_int.shape[0]
    n_allzero = int((counts_int.sum(axis=1) == 0).sum())

    keep = counts_int.sum(axis=1) >= 10
    counts_f = counts_int[keep]

    print(f"genes: {n_before} in file -> {counts_f.shape[0]} retained "
          f"({100 * counts_f.shape[0] / n_before:.1f}%)")
    print(f"all-zero across all samples: {n_allzero} "
          f"({100 * n_allzero / n_before:.1f}%)")

    # A parsing/filtering regression would show up immediately as a collapsed gene count.
    assert counts_f.shape[0] > 15000, "far fewer genes retained than expected"
    return counts_f


def run_deseq(counts_f: pd.DataFrame, coldata: pd.DataFrame):
    """
    Build a DeseqDataSet and fit it. ``counts`` for PyDESeq2 is samples x genes (AnnData
    orientation), so the gene x sample matrix is transposed here.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference

    counts_sxg = counts_f.T  # samples x genes
    metadata = coldata.set_index("ID")[["cellLine", "treatment"]].loc[counts_sxg.index]
    metadata["cellLine"] = metadata["cellLine"].astype(str)
    # Set the reference level as the FIRST category: formulaic uses the first category of a
    # categorical as the base, so this makes Untreated the reference and yields a clean
    # 'treatment[T.Dexamethasone]' coefficient -- positive LFC = induced by dexamethasone.
    # (PyDESeq2's ref_level argument is deprecated in 0.5.x; encoding it here is the
    # supported route, and the coefficient-name assertion below fails loudly if it slips.)
    metadata["treatment"] = pd.Categorical(
        metadata["treatment"], categories=["Untreated", "Dexamethasone"]
    )
    assert list(metadata.index) == list(counts_sxg.index), "metadata/count sample mismatch"

    dds = DeseqDataSet(
        counts=counts_sxg,
        metadata=metadata,
        design="~cellLine + treatment",
        refit_cooks=True,
        inference=DefaultInference(),
    )
    dds.deseq2()
    return dds


def extract_results(dds, symbol_of):
    """Run the Wald test for the treatment contrast and merge shrunken fold changes."""
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference

    rule("3. Differential expression: Dexamethasone vs Untreated")

    stat_res = DeseqStats(dds, contrast=COEF_CONTRAST, inference=DefaultInference())
    stat_res.summary()
    res_unshrunk = stat_res.results_df.copy()

    # Determine the coefficient name from the fitted model rather than hard-coding it, then
    # assert it exists — the same "can't extract the wrong contrast" guard as the R version.
    lfc_cols = list(dds.varm["LFC"].columns)
    coeff = next((c for c in lfc_cols if "treatment" in c and "Dexamethasone" in c), None)
    assert coeff is not None, f"no Dexamethasone treatment coefficient in {lfc_cols}"
    print(f"shrinking coefficient: {coeff}")

    stat_res.lfc_shrink(coeff=coeff)
    res_shrunk = stat_res.results_df.copy()

    # Carry BOTH: unshrunken supplies p-values/padj and raw effect sizes; shrunken supplies
    # the effect sizes used for ranking and plotting. Keeping them side by side avoids the
    # common error of mixing the two.
    res = res_unshrunk.copy()
    res["log2FoldChange_shrunk"] = res_shrunk["log2FoldChange"]
    res["lfcSE_shrunk"] = res_shrunk["lfcSE"]
    res = res.reset_index().rename(columns={"index": "gene_id"})
    if "gene_id" not in res.columns:  # PyDESeq2 index name varies by version
        res = res.rename(columns={res.columns[0]: "gene_id"})

    res.insert(1, "symbol", [symbol_of(g) for g in res["gene_id"]])
    res = res[[
        "gene_id", "symbol", "baseMean",
        "log2FoldChange", "lfcSE",
        "log2FoldChange_shrunk", "lfcSE_shrunk",
        "stat", "pvalue", "padj",
    ]]
    res = res.sort_values(
        ["padj", "log2FoldChange_shrunk"],
        key=lambda s: -s.abs() if s.name == "log2FoldChange_shrunk" else s,
        na_position="last",
    ).reset_index(drop=True)

    sig = res["padj"] < PADJ_CUTOFF
    n_sig = int(sig.sum())
    n_up = int((sig & (res["log2FoldChange_shrunk"] > 0)).sum())
    n_down = int((sig & (res["log2FoldChange_shrunk"] < 0)).sum())
    n_strong = int((sig & (res["log2FoldChange_shrunk"].abs() > LFC_CUTOFF)).sum())

    print(f"genes tested (padj not NA): {int(res['padj'].notna().sum())}")
    print(f"padj < {PADJ_CUTOFF}: {n_sig}  ({n_up} up, {n_down} down)")
    print(f"padj < {PADJ_CUTOFF} & |LFC| > {LFC_CUTOFF}: {n_strong}")

    return res, n_sig


def variance_stabilise(dds):
    """VST for QC/clustering, with a log2(normalized+1) fallback across PyDESeq2 versions."""
    genes = list(dds.var_names)
    samples = list(dds.obs_names)
    try:
        try:
            dds.vst(use_design=True)  # design-aware, mirrors R's blind = FALSE
        except TypeError:
            dds.vst()
        mat = np.asarray(dds.layers["vst_counts"])
        name = "VST"
    except Exception:
        normed = np.asarray(dds.layers["normed_counts"])
        mat = np.log2(normed + 1)
        name = "log2(normalized+1)"
    return pd.DataFrame(mat, index=samples, columns=genes), name


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------


def plot_library_sizes(dds, coldata, figdir):
    totals = np.asarray(dds.X).sum(axis=1)
    df = pd.DataFrame({
        "sample": list(dds.obs_names),
        "total": totals / 1e6,
        "treatment": list(dds.obs["treatment"]),
    }).sort_values("total")
    colours = {"Untreated": "#66c2a5", "Dexamethasone": "#fc8d62"}

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.barh(df["sample"], df["total"],
            color=[colours.get(t, "grey") for t in df["treatment"]])
    ax.set_xlabel("Total counts (millions)")
    ax.set_title("Library sizes")
    ax.legend(handles=[Line2D([0], [0], marker="s", linestyle="", color=c, label=t)
                       for t, c in colours.items()], title="Treatment")
    savefig(fig, figdir, "fig-library-sizes")
    print(f"  library sizes: {df['total'].min():.1f}M - {df['total'].max():.1f}M "
          f"({df['total'].max() / df['total'].min():.2f}-fold spread)")


def plot_pca(vst_df, dds, figdir):
    scores, pct = pca_2d(vst_df.values)
    meta = dds.obs
    cell_lines = list(meta["cellLine"])
    treatments = list(meta["treatment"])

    uniq_cl = sorted(set(cell_lines))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    marker_of = {cl: markers[i % len(markers)] for i, cl in enumerate(uniq_cl)}
    colours = {"Untreated": "#377eb8", "Dexamethasone": "#e41a1c"}

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    # Dashed connectors between the paired samples of each cell line.
    for cl in uniq_cl:
        idx = [i for i, c in enumerate(cell_lines) if c == cl]
        ax.plot(scores[idx, 0], scores[idx, 1], color="grey", ls="--", lw=1, zorder=1)
    for i in range(scores.shape[0]):
        ax.scatter(scores[i, 0], scores[i, 1], s=90, zorder=2,
                   color=colours.get(treatments[i], "grey"),
                   marker=marker_of[cell_lines[i]], edgecolor="black", linewidth=0.4)
    ax.set_xlabel(f"PC1: {round(100 * pct[0])}% variance")
    ax.set_ylabel(f"PC2: {round(100 * pct[1])}% variance")
    ax.set_title("Samples cluster by cell line; treatment shifts within pairs")

    handles = [Line2D([0], [0], marker="o", linestyle="", color=c, label=t)
               for t, c in colours.items()]
    handles += [Line2D([0], [0], marker=marker_of[cl], linestyle="", color="grey", label=cl)
                for cl in uniq_cl]
    ax.legend(handles=handles, fontsize=8, loc="best")
    savefig(fig, figdir, "fig-pca")


def plot_sample_distances(vst_df, dds, figdir):
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import linkage, leaves_list

    d = squareform(pdist(vst_df.values))
    order = leaves_list(linkage(pdist(vst_df.values), method="average"))
    labels = [f"{cl} / {tr}" for cl, tr in zip(dds.obs["cellLine"], dds.obs["treatment"])]
    labels = np.array(labels)

    d = d[np.ix_(order, order)]
    lab = labels[order]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(d, cmap="Blues_r")
    ax.set_xticks(range(len(lab)))
    ax.set_yticks(range(len(lab)))
    ax.set_xticklabels(lab, rotation=90, fontsize=8)
    ax.set_yticklabels(lab, fontsize=8)
    ax.set_title("Sample distances (VST)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    savefig(fig, figdir, "fig-sample-distances")


def plot_dispersions(dds, figdir):
    try:
        # PyDESeq2 stores per-gene dispersions and mean counts in dds.var (0.5.x).
        means = np.asarray(dds.var["_normed_means"])
        genewise = np.asarray(dds.var["genewise_dispersions"])
        fitted = np.asarray(dds.var["fitted_dispersions"])
        final = np.asarray(dds.var["dispersions"])
        ok = (means > 0) & np.isfinite(genewise)

        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.scatter(means[ok], genewise[ok], s=4, color="black", alpha=0.3, label="genewise")
        ax.scatter(means[ok], final[ok], s=4, color="#377eb8", alpha=0.3, label="final (MAP)")
        o = np.argsort(means[ok])
        ax.plot(means[ok][o], fitted[ok][o], color="#e41a1c", lw=1.5, label="fitted trend")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("mean of normalized counts")
        ax.set_ylabel("dispersion")
        ax.set_title("Dispersion estimates")
        ax.legend()
        savefig(fig, figdir, "fig-dispersion")
    except Exception as exc:
        print(f"  dispersion plot skipped ({type(exc).__name__}: {exc})")


def plot_ma(res, figdir):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    d = res[res["baseMean"] > 0]
    sig = d["padj"] < PADJ_CUTOFF
    ax.scatter(d["baseMean"][~sig], d["log2FoldChange_shrunk"][~sig],
               s=5, color="grey", alpha=0.4, label="ns")
    ax.scatter(d["baseMean"][sig], d["log2FoldChange_shrunk"][sig],
               s=6, color="#e6550d", alpha=0.6, label=f"padj<{PADJ_CUTOFF}")
    ax.set_xscale("log")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(-5, 5)
    ax.set_xlabel("mean of normalized counts")
    ax.set_ylabel("log2 fold change (shrunken)")
    ax.set_title("MA plot")
    ax.legend()
    savefig(fig, figdir, "fig-ma")


def plot_volcano(res, label_of, figdir):
    d = res[res["padj"].notna()].copy()
    strong = (d["padj"] < PADJ_CUTOFF) & (d["log2FoldChange_shrunk"].abs() > LFC_CUTOFF)
    mid = (d["padj"] < PADJ_CUTOFF) & ~strong
    ns = d["padj"] >= PADJ_CUTOFF

    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(d["log2FoldChange_shrunk"][ns], -np.log10(d["padj"][ns]),
               s=6, color="grey", alpha=0.5, label="not significant")
    ax.scatter(d["log2FoldChange_shrunk"][mid], -np.log10(d["padj"][mid]),
               s=8, color="#3182bd", alpha=0.6, label=f"padj<{PADJ_CUTOFF}")
    ax.scatter(d["log2FoldChange_shrunk"][strong], -np.log10(d["padj"][strong]),
               s=8, color="#e6550d", alpha=0.7,
               label=f"padj<{PADJ_CUTOFF} & |LFC|>{LFC_CUTOFF}")
    ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", color="grey", lw=0.8)
    for x in (-LFC_CUTOFF, LFC_CUTOFF):
        ax.axvline(x, ls="--", color="grey", lw=0.8)

    for _, row in d.nsmallest(15, "padj").iterrows():
        ax.annotate(label_of(row["gene_id"]),
                    (row["log2FoldChange_shrunk"], -np.log10(row["padj"])),
                    fontsize=7, ha="center", va="bottom")
    ax.set_xlabel("log2 fold change (Dexamethasone vs Untreated, shrunken)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Dexamethasone response")
    ax.legend(fontsize=8, loc="upper center")
    savefig(fig, figdir, "fig-volcano")


def plot_pvalue_hist(res, figdir):
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.hist(res["pvalue"].dropna(), bins=50, color="#3182bd", edgecolor="white")
    ax.set_xlabel("raw p-value")
    ax.set_ylabel("genes")
    ax.set_title("p-value distribution")
    savefig(fig, figdir, "fig-pvalue-hist")


def plot_top_gene_heatmap(res, vst_df, dds, label_of, figdir):
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    top = res.loc[res["padj"].notna(), "gene_id"].head(40).tolist()
    mat = vst_df[top].T.values  # genes x samples
    mat = mat - mat.mean(axis=1, keepdims=True)  # centre per gene

    col_order = leaves_list(linkage(pdist(mat.T), method="average"))
    row_order = leaves_list(linkage(pdist(mat), method="average"))
    mat = mat[np.ix_(row_order, col_order)]

    gene_labels = [label_of(top[i]) for i in row_order]
    samples = list(dds.obs_names)
    col_labels = [f"{cl} / {tr}" for cl, tr in zip(dds.obs["cellLine"], dds.obs["treatment"])]
    col_labels = [col_labels[i] for i in col_order]

    vmax = np.abs(mat).max()
    fig, ax = plt.subplots(figsize=(7, 9))
    im = ax.imshow(mat, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(gene_labels)))
    ax.set_yticklabels(gene_labels, fontsize=7)
    ax.set_title("Top 40 DE genes (centred VST)")
    fig.colorbar(im, ax=ax, shrink=0.5)
    savefig(fig, figdir, "fig-top-genes-heatmap")


def plot_markers(dds, res, have_symbols, figdir):
    normed = pd.DataFrame(np.asarray(dds.layers["normed_counts"]),
                          index=dds.obs_names, columns=dds.var_names)
    cell_lines = list(dds.obs["cellLine"])
    treatments = list(dds.obs["treatment"])
    t_levels = ["Untreated", "Dexamethasone"]

    available = {sym: gid for sym, gid in MARKER_GENES.items() if gid in normed.columns}
    missing = [sym for sym in MARKER_GENES if sym not in available]
    if missing:
        print(f"  markers skipped (not in filtered set): {', '.join(missing)}")
    if not available:
        print("  no marker genes available to plot")
        return

    uniq_cl = sorted(set(cell_lines))
    cmap = plt.get_cmap("tab10")
    colour_of = {cl: cmap(i) for i, cl in enumerate(uniq_cl)}

    n = len(available)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), squeeze=False)

    for k, (sym, gid) in enumerate(available.items()):
        ax = axes[k // ncol][k % ncol]
        counts = normed[gid]
        for cl in uniq_cl:
            xs, ys = [], []
            for lvl_i, lvl in enumerate(t_levels):
                vals = [counts[s] for s, c, t in zip(normed.index, cell_lines, treatments)
                        if c == cl and t == lvl]
                for v in vals:
                    xs.append(lvl_i)
                    ys.append(v)
            order = np.argsort(xs)
            ax.plot(np.array(xs)[order], np.array(ys)[order],
                    color="grey", ls="--", lw=1, zorder=1)
            ax.scatter(xs, ys, color=colour_of[cl], s=40, zorder=2, label=cl)
        ax.set_yscale("log")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(t_levels, rotation=15)
        ax.set_title(sym)
        ax.set_ylabel("normalized count")

    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    handles = [Line2D([0], [0], marker="o", linestyle="", color=colour_of[cl], label=cl)
               for cl in uniq_cl]
    fig.legend(handles=handles, title="Cell line", loc="lower center",
               ncol=len(uniq_cl), fontsize=8)
    fig.suptitle("Known dexamethasone responders")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    savefig(fig, figdir, "fig-marker-genes")

    marker_stats = res[res["gene_id"].isin(available.values())]
    print("  marker gene statistics (positive log2FC = induced by dexamethasone):")
    for _, row in marker_stats.iterrows():
        sym = next(s for s, g in available.items() if g == row["gene_id"])
        print(f"    {sym:>9} {row['gene_id']}  "
              f"log2FC={row['log2FoldChange_shrunk']:+.2f}  padj={row['padj']:.2e}")


# --------------------------------------------------------------------------------------
# Outputs and diagnostics
# --------------------------------------------------------------------------------------


def write_tables(res, dds, vst_df, symbol_of, outdir):
    rule("4. Export results")

    sig = res[res["padj"].notna() & (res["padj"] < PADJ_CUTOFF)]

    normed = pd.DataFrame(np.asarray(dds.layers["normed_counts"]).T,
                          index=dds.var_names, columns=dds.obs_names)
    normed.insert(0, "symbol", [symbol_of(g) for g in normed.index])
    normed.insert(0, "gene_id", normed.index)

    vst_out = vst_df.T.copy()
    vst_out.insert(0, "symbol", [symbol_of(g) for g in vst_out.index])
    vst_out.insert(0, "gene_id", vst_out.index)

    sf = pd.DataFrame({
        "sample": list(dds.obs_names),
        "size_factor": np.asarray(dds.obs["size_factors"]),
        "library_size": np.asarray(dds.X).sum(axis=1),
        "cellLine": list(dds.obs["cellLine"]),
        "treatment": list(dds.obs["treatment"]),
    })

    tables = {
        "de_results_all.tsv": res,
        "de_results_significant.tsv": sig,
        "normalized_counts.tsv": normed,
        "vst_counts.tsv": vst_out,
        "size_factors.tsv": sf,
    }
    for name, df in tables.items():
        df.to_csv(outdir / name, sep="\t", index=False, na_rep="NA")
        print(f"  results/{name:<32} {len(df):6d} rows")


def paired_vs_unpaired(counts_f, coldata, res, n_sig_paired, figdir):
    """Refit ~treatment alone to quantify what discarding the pairing costs."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference

    rule("5. Diagnostic — paired vs unpaired")

    counts_sxg = counts_f.T
    metadata = coldata.set_index("ID")[["cellLine", "treatment"]].loc[counts_sxg.index]
    metadata["treatment"] = pd.Categorical(
        metadata["treatment"], categories=["Untreated", "Dexamethasone"]
    )
    inference = DefaultInference()

    dds_u = DeseqDataSet(counts=counts_sxg, metadata=metadata,
                         design="~treatment", refit_cooks=True, inference=inference)
    dds_u.deseq2()

    stat_u = DeseqStats(dds_u, contrast=COEF_CONTRAST, inference=inference)
    stat_u.summary()
    res_u = stat_u.results_df

    n_sig_unpaired = int((res_u["padj"] < PADJ_CUTOFF).sum())
    print(f"  ~ cellLine + treatment (paired):  {n_sig_paired} genes at padj<{PADJ_CUTOFF}")
    print(f"  ~ treatment (unpaired):           {n_sig_unpaired} genes at padj<{PADJ_CUTOFF}")
    print(f"  unpaired recovers {100 * n_sig_unpaired / max(n_sig_paired, 1):.0f}% of the paired hits")

    paired_padj = res.set_index("gene_id")["padj"]
    cmp = pd.DataFrame({
        "unpaired": -np.log10(res_u["padj"]),
        "paired": -np.log10(paired_padj.reindex(res_u.index)),
    }).replace([np.inf, -np.inf], np.nan).dropna()

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(cmp["unpaired"], cmp["paired"], s=5, alpha=0.25)
    lim = max(cmp["unpaired"].max(), cmp["paired"].max())
    ax.plot([0, lim], [0, lim], color="#e6550d")
    ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", color="grey", lw=0.8)
    ax.axvline(-np.log10(PADJ_CUTOFF), ls="--", color="grey", lw=0.8)
    ax.set_xlabel("-log10 padj, unpaired (~ treatment)")
    ax.set_ylabel("-log10 padj, paired (~ cellLine + treatment)")
    ax.set_title("Blocking on cell line recovers evidence the unpaired model loses")
    savefig(fig, figdir, "fig-paired-vs-unpaired")


# --------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Airway RNA-seq DE analysis (PyDESeq2).")
    repo = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=repo / "data")
    parser.add_argument("--results-dir", type=Path, default=repo / "results")
    parser.add_argument("--figures-dir", type=Path, default=repo / "figures")
    parser.add_argument("--ensembl-release", type=int, default=110,
                        help="Ensembl release for pyensembl symbol lookup (optional).")
    args = parser.parse_args(argv)

    np.random.seed(RANDOM_SEED)
    args.results_dir.mkdir(exist_ok=True)
    args.figures_dir.mkdir(exist_ok=True)

    try:
        import pydeseq2  # noqa: F401
    except ImportError:
        print("PyDESeq2 is required. Install with:\n"
              "  pip install -r setup/requirements.txt", file=sys.stderr)
        return 1

    # 1. Load + validate
    coldata, counts_gxs = load_inputs(args.data_dir)
    show_design(coldata)

    # 2. Prepare counts
    counts_f = prepare_counts(counts_gxs)

    # Symbols (optional, offline)
    rule("Gene symbols (optional)")
    symbol_map, have_symbols = build_symbol_map(list(counts_f.index), args.ensembl_release)
    if have_symbols:
        mapped = sum(1 for v in symbol_map.values() if v)
        print(f"  symbols mapped: {mapped}/{len(symbol_map)} "
              f"({100 * mapped / len(symbol_map):.1f}%)")
    label_of, symbol_of = make_labeller(symbol_map)

    # 3. Fit + results
    dds = run_deseq(counts_f, coldata)
    res, n_sig = extract_results(dds, symbol_of)

    top = res.head(20).copy()
    top["gene"] = [label_of(g) for g in top["gene_id"]]
    print("\ntop 20 genes by padj:")
    print(top[["gene", "gene_id", "baseMean", "log2FoldChange_shrunk", "padj"]]
          .to_string(index=False))

    # QC transform + figures
    rule("Figures")
    vst_df, transform_name = variance_stabilise(dds)
    print(f"  QC transform: {transform_name}")
    plot_library_sizes(dds, coldata, args.figures_dir)
    plot_pca(vst_df, dds, args.figures_dir)
    plot_sample_distances(vst_df, dds, args.figures_dir)
    plot_dispersions(dds, args.figures_dir)
    plot_ma(res, args.figures_dir)
    plot_volcano(res, label_of, args.figures_dir)
    plot_pvalue_hist(res, args.figures_dir)
    plot_top_gene_heatmap(res, vst_df, dds, label_of, args.figures_dir)
    plot_markers(dds, res, have_symbols, args.figures_dir)
    print(f"  figures written to {args.figures_dir}")

    # 4. Tables
    write_tables(res, dds, vst_df, symbol_of, args.results_dir)

    # 5. Diagnostic
    paired_vs_unpaired(counts_f, coldata, res, n_sig, args.figures_dir)

    rule("Done")
    print(f"results: {args.results_dir}")
    print(f"figures: {args.figures_dir}")
    print(f"pandas {pd.__version__}, numpy {np.__version__}, "
          f"pydeseq2 {pydeseq2.__version__ if hasattr(pydeseq2, '__version__') else '?'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
