#!/usr/bin/env python3
"""Build a gene-of-interest isoform table from IsoQuant grouped counts + SQANTI3
classification. Joins per-platform counts with the isoform's structural category
and associated gene, then keeps rows for the target gene (default APP).

SQANTI3 reports `associated_gene` as gene IDs (e.g. ENSG00000142192.22; fusions join
several with "_"), so the symbol is resolved to its gene ID(s) through the reference
GTF and matched exactly (version ignored).
"""
import argparse
import gzip
import re

import pandas as pd


def load_classification(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", dtype=str)
        frames.append(df)
    cls = pd.concat(frames, ignore_index=True)
    # SQANTI3 columns: isoform, associated_gene, structural_category, ...
    keep = [c for c in ["isoform", "associated_gene", "structural_category"] if c in cls.columns]
    return cls[keep]


def load_counts(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t")
        # first column is the feature/transcript id
        df = df.rename(columns={df.columns[0]: "isoform"})
        frames.append(df.set_index("isoform"))
    counts = pd.concat(frames, axis=1)
    return counts.reset_index()


def gene_ids_for(symbol, gtf):
    """Unversioned gene IDs whose gene_name equals `symbol` in the GTF."""
    ids = set()
    opener = gzip.open if gtf.endswith(".gz") else open
    pat = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "gene":
                continue
            m = pat.search(f[8])
            if m and m.group(2) == symbol:
                ids.add(m.group(1).split(".")[0])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="APP")
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--classification", nargs="+", required=True)
    ap.add_argument("--gtf", required=True, help="reference GTF, to resolve the gene symbol to gene IDs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cls = load_classification(args.classification)
    counts = load_counts(args.counts)
    ids = gene_ids_for(args.gene, args.gtf)
    if not ids:
        raise SystemExit(f"[app_isoform_table] gene {args.gene} not found in {args.gtf}")

    merged = cls.merge(counts, on="isoform", how="inner")
    assoc = merged["associated_gene"].fillna("").str.split("_")
    hit = assoc.apply(lambda genes: any(g.split(".")[0] in ids for g in genes))
    sub = merged[hit]
    sub.to_csv(args.out, sep="\t", index=False)
    print(f"[app_isoform_table] gene={args.gene} ({','.join(sorted(ids))}) rows={len(sub)} -> {args.out}")


if __name__ == "__main__":
    main()
