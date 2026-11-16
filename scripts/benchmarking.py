#!/usr/bin/env python3
"""Collect Snakemake benchmark files into one table (runtime, memory, I/O per rule).

Usage: benchmarking.py OUTDIR > benchmarking.tsv
Reads OUTDIR/<group>/benchmarks/<group>_<rule>.tsv for every group directory.
"""
import glob
import os
import sys

import pandas as pd

COLS = ["sample", "rule", "s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in", "io_out", "mean_load"]


def hms(seconds: float) -> str:
    m, s = divmod(int(round(seconds)), 60)
    h, m = divmod(m, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


def collect(outdir: str) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(os.path.join(outdir, "*", "benchmarks", "*.tsv"))):
        group = os.path.basename(os.path.dirname(os.path.dirname(f)))
        table = pd.read_csv(f, sep="\t")
        for c in COLS[2:]:
            if c in table.columns:
                table[c] = pd.to_numeric(table[c], errors="coerce")
        table["sample"] = group
        table["rule"] = os.path.splitext(os.path.basename(f))[0].replace(group + "_", "", 1)
        rows.append(table[[c for c in COLS if c in table.columns]])
    if not rows:
        return pd.DataFrame(columns=COLS)
    df = pd.concat(rows)
    agg = df.groupby(["sample", "rule"]).agg("max")
    agg["s"] = agg["s"].apply(hms)
    return agg


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    collect(sys.argv[1]).to_csv(sys.stdout, sep="\t")
