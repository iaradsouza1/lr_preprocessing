#!/usr/bin/env python3
"""Build a gene-of-interest isoform table from IsoQuant grouped counts + SQANTI3
classification. Joins per-platform counts with the isoform's structural category
and associated gene, then keeps rows for the target gene (default APP).
"""
import argparse

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", default="APP")
    ap.add_argument("--counts", nargs="+", required=True)
    ap.add_argument("--classification", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cls = load_classification(args.classification)
    counts = load_counts(args.counts)

    merged = cls.merge(counts, on="isoform", how="inner")
    if "associated_gene" in merged.columns:
        sub = merged[merged["associated_gene"].astype(str).str.contains(args.gene, case=False, na=False)]
    else:
        sub = merged
    sub.to_csv(args.out, sep="\t", index=False)
    print(f"[app_isoform_table] gene={args.gene} rows={len(sub)} -> {args.out}")


if __name__ == "__main__":
    main()
