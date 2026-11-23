"""Merge per-sample IsoQuant requant tables (isoquant_requant_unfiltered.py) into one matrix.

Transcripts confirmed in no sample are dropped, as IsoQuant does in a joint run; TPM is then
count / column total * 1e6 (IsoQuant's default "simple" normalization).
Usage: merge_requant.py OUT_COUNTS OUT_TPM --counts S1.tsv ... --confirmed S1.txt ...
"""
import argparse

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("out_counts")
p.add_argument("out_tpm")
p.add_argument("--counts", nargs="+", required=True)
p.add_argument("--confirmed", nargs="+", required=True)
a = p.parse_args()

confirmed = set()
for f in a.confirmed:
    confirmed.update(line.strip() for line in open(f) if line.strip())
tabs = [pd.read_csv(f, sep="\t", index_col=0) for f in a.counts]
counts = pd.concat(tabs, axis=1).fillna(0.0)
counts.index.name = tabs[0].index.name
# a joint run writes only confirmed transcripts (the others are zeroed and dropped)
counts = counts.loc[counts.index.isin(confirmed)].sort_index()
counts.to_csv(a.out_counts, sep="\t", float_format="%.2f")
(counts / counts.sum() * 1e6).to_csv(a.out_tpm, sep="\t", float_format="%.6f")
print(f"{len(counts)} confirmed transcripts, {len(a.counts)} samples")
