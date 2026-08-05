#!/usr/bin/env python3
"""
Pathway enrichment via gseapy's Enrichr client — the working fallback for
clawbio's pathway-enricher (whose addList POST is form-encoded and 400s on
real gene lists; Enrichr requires multipart, which gseapy uses correctly).

Usage:  python enrich_gseapy.py <gene_list.txt> <output_dir>
Writes: <output_dir>/gseapy_enrichment.csv   (all libraries, long format)
"""
import os
import sys

import gseapy as gp

LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: enrich_gseapy.py <gene_list.txt> <output_dir>", file=sys.stderr)
        return 2
    inp, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    genes = [ln.strip() for ln in open(inp) if ln.strip()]
    if len(genes) < 5:
        print(f"too few genes ({len(genes)}) for enrichment", file=sys.stderr)
        return 1

    enr = gp.enrichr(gene_list=genes, gene_sets=LIBRARIES, outdir=None, no_plot=True)
    res = enr.results.rename(
        columns={
            "Gene_set": "library",
            "P-value": "pvalue",
            "Adjusted P-value": "adj_pvalue",
            "Combined Score": "combined_score",
        }
    )
    out_csv = os.path.join(outdir, "gseapy_enrichment.csv")
    res.to_csv(out_csv, index=False)
    print(f"gseapy: {len(res)} rows across {res['library'].nunique()} libraries -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
