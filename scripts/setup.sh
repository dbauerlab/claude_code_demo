#!/usr/bin/env bash
#
# setup.sh — one-shot toolchain install for the airway RNA-seq DE workflow.
#
# Installs (globally, WITHOUT sudo/admin password):
#   - Quarto            user-space tarball  (renders airway_de.qmd)
#   - R + R packages    Homebrew formula    (rendering, plotting, tables, Ensembl->symbol)
#   - Python 3.12       Homebrew formula    (clawbio needs 3.11-3.12; the machine's 3.14/3.9 won't do)
#   - Python packages   pip                 (pydeseq2 + the clawbio rnaseq-de / pathway-enricher stack)
#
# The differential-expression ENGINE is clawbio (PyDESeq2), not R/DESeq2 — so
# DESeq2/Bioconductor DE packages are deliberately NOT installed. R only renders/visualizes.
#
# Sudo-free by design: R & Python use Homebrew *formulae* (installed under the
# Homebrew prefix, no admin prompt), and Quarto installs from its release tarball
# into ~/.local. The `--cask` installers for R/Quarto need an admin password and
# are avoided here; if you prefer them, install by hand from quarto.org / CRAN.
#
# Usage:  bash scripts/setup.sh
#
set -euo pipefail
export NONINTERACTIVE=1 HOMEBREW_NO_AUTO_UPDATE=1

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m[warn] %s\033[0m\n' "$*"; }

# --- Homebrew on PATH (Apple Silicon default) --------------------------------
if [ -x /opt/homebrew/bin/brew ]; then
  eval "$(/opt/homebrew/bin/brew shellenv)"
elif ! command -v brew >/dev/null 2>&1; then
  warn "Homebrew not found. Install it from https://brew.sh, then re-run."
  exit 1
fi
BREW_BIN="$(brew --prefix)/bin"

# --- Quarto (user-space tarball, no sudo) ------------------------------------
if command -v quarto >/dev/null 2>&1; then
  say "Quarto already installed ($(quarto --version))"
else
  say "Installing Quarto (user-space tarball)"
  QVER="$(curl -fsSLI -o /dev/null -w '%{url_effective}' \
          https://github.com/quarto-dev/quarto-cli/releases/latest | sed 's#.*/tag/v##')"
  [ -n "$QVER" ] || { warn "Could not resolve latest Quarto version"; exit 1; }
  TARBALL="quarto-${QVER}-macos.tar.gz"
  curl -fsSL -o "/tmp/${TARBALL}" \
    "https://github.com/quarto-dev/quarto-cli/releases/download/v${QVER}/${TARBALL}"
  DEST="$HOME/.local/quarto-${QVER}"
  mkdir -p "$DEST"
  tar -xzf "/tmp/${TARBALL}" -C "$DEST" --strip-components=1
  ln -sf "${DEST}/bin/quarto" "${BREW_BIN}/quarto"   # BREW_BIN is user-writable & on PATH
  say "Quarto $(quarto --version) installed"
fi

# --- R (Homebrew formula, no sudo) -------------------------------------------
if command -v Rscript >/dev/null 2>&1; then
  say "R already installed ($(Rscript -e 'cat(R.version.string)'))"
else
  say "Installing R (brew formula)"
  brew install r
fi

# --- R packages (rendering + viz + Ensembl->symbol) --------------------------
say "Installing R packages (CRAN + Bioconductor org.Hs.eg.db) — this can take a while"
Rscript - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1, parallel::detectCores()))
cran <- c("readr","dplyr","tidyr","tibble","ggplot2","ggrepel",
          "pheatmap","RColorBrewer","jsonlite","DT","knitr","rmarkdown")
need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

cat("R packages OK\n")
RS

# --- Python 3.12 + clawbio deps ----------------------------------------------
if command -v python3.12 >/dev/null 2>&1; then
  say "python3.12 already installed ($(python3.12 --version))"
else
  say "Installing python@3.12 (brew formula)"
  brew install python@3.12
fi

say "Installing Python packages (pydeseq2 + clawbio stack) into python3.12"
# NB: do NOT `pip install --upgrade pip` here — Homebrew's pip ships without a
# RECORD file, so a self-upgrade fails with 'uninstall-no-record-file'. Not needed.
python3.12 -m pip install --break-system-packages \
  "pydeseq2" "pandas" "numpy" "scipy" "scikit-learn" "matplotlib" \
  "gseapy" "requests" \
  "opentelemetry-sdk>=1.20,<2" "rocrate>=0.15.0" "biopython>=1.82" "pyyaml>=6.0"
  # ^ the last line are clawbio-package runtime deps its skill scripts import at load time

# --- Final check -------------------------------------------------------------
say "Verifying toolchain"
quarto --version
Rscript -e 'suppressMessages(library(org.Hs.eg.db)); cat("R + org.Hs.eg.db OK\n")'
python3.12 -c "import pydeseq2, pandas, sklearn, scipy, matplotlib; print('Python DE stack OK')"
python3.12 -c "import gseapy; print('gseapy (enrichment) OK')" || warn "gseapy import failed (enrichment degrades gracefully)"

say "Done. Render the analysis with:  quarto render airway_de.qmd"
