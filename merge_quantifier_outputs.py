#!/usr/bin/env python3
"""
Merge per-sample miRDeep2 quantifier outputs into one raw count matrix.

Run this from the folder that contains the per-sample directories (XS1, XS2, ...).

quantifier.pl was run with -W, so a read mapping to N precursors contributes
1/N to each. Summing a mature miRNA's precursor rows is therefore correct.
Without -W each row would carry the full count and summing would inflate it.
"""

import os
import sys

import pandas as pd

ROOT = "."                              # folder containing the per-sample directories
SAMPLES = [f"XS{i}" for i in range(1, 13)]
OUTPUT = "mirna_count_matrix.tsv"


def load(sample):
    """One sample -> a Series of counts per mature miRNA."""
    path = os.path.join(ROOT, sample, f"miRNAs_expressed_all_samples_{sample}.csv")
    if not os.path.exists(path):
        sys.exit(f"ERROR: not found: {path}\n"
                 f"Run this from the folder containing {SAMPLES[0]}/, or set ROOT.")
    print(f"{sample} -> {path}")
    return pd.read_csv(path, sep="\t").groupby("#miRNA")["read_count"].sum()


df = pd.DataFrame({s: load(s) for s in SAMPLES})

# quantifier.pl emits a row for every database entry, so gaps should not happen.
if df.isna().any().any():
    print(f"WARNING: {int(df.isna().sum().sum())} miRNA/sample cells absent -> 0. "
          "Did every sample use the same database?")

df = df.fillna(0).round().astype(int).sort_index()   # -W gives fractional counts
df.to_csv(OUTPUT, sep="\t", index_label="miRNA")

print(f"\nWrote {OUTPUT}: {df.shape[0]} miRNAs x {df.shape[1]} samples")
print("Library sizes:")
print(df.sum().to_string())
