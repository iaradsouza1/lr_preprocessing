#!/usr/bin/env python3
"""QC gate 2: mean minimap2 de:f (per-base divergence) identity from a BAM.

identity = 1 - de  (de is minimap2's gap-compressed per-base divergence tag).
Reports mean identity over primary alignments and compares to a threshold.
"""
import argparse
import sys

import pysam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--gate", choices=["warn", "fail"], default="warn")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    total = 0.0
    n = 0
    missing = 0
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            try:
                de = read.get_tag("de")
            except KeyError:
                missing += 1
                continue
            total += 1.0 - float(de)
            n += 1

    mean_identity = (total / n) if n else float("nan")
    passed = n > 0 and mean_identity >= args.threshold

    with open(args.out, "w") as fh:
        fh.write(f"bam\t{args.bam}\n")
        fh.write(f"n_primary\t{n}\n")
        fh.write(f"n_missing_de_tag\t{missing}\n")
        fh.write(f"mean_identity\t{mean_identity:.5f}\n")
        fh.write(f"threshold\t{args.threshold}\n")
        fh.write(f"pass\t{passed}\n")

    msg = f"[identity_check] {args.bam}: mean_identity={mean_identity:.5f} threshold={args.threshold} pass={passed}"
    print(msg, file=sys.stderr)
    if not passed and args.gate == "fail":
        sys.exit(f"IDENTITY GATE FAILED: {msg}")


if __name__ == "__main__":
    main()
