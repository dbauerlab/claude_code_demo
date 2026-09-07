#!/usr/bin/env bash
set -euo pipefail
ROOT="$(git rev-parse --show-toplevel)"
DEST="$ROOT/analysis/work/branches"
mkdir -p "$DEST"
BRANCHES="Catarina Olivia santosnuno okp-23 becky tiph barretj pateld1 davidlvb marga"
for b in $BRANCHES; do
  d="$DEST/$b"
  rm -rf "$d"; mkdir -p "$d"
  git archive "origin/demo/$b" | tar -x -C "$d"
  echo "materialized $b -> $(ls "$d" | tr '\n' ' ')"
done
