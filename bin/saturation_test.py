#!/usr/bin/env python3
"""
Alignment QC statistics and Xia's substitution saturation test.

Computes:
  - Alignment statistics: N sequences, length, gap%, informative sites
  - Xia's Iss (index of substitution saturation) per codon position
  - Saturation scatter plot (transitions vs transversions)

Usage:
  saturation_test.py --input gene.trimmed.fasta --gene blaCTX-M \
                     --out-stats stats.csv --out-saturation sat.csv --out-plot sat.png
"""

import argparse
import csv
import sys
from collections import Counter
from itertools import combinations

import numpy as np
from Bio import SeqIO

# Try importing matplotlib; skip plot if not available
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# Transition pairs (purine ↔ purine, pyrimidine ↔ pyrimidine)
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_args():
    parser = argparse.ArgumentParser(description="Alignment QC + saturation test")
    parser.add_argument("--input", required=True, help="Input alignment FASTA")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--out-stats", required=True, help="Output alignment stats CSV")
    parser.add_argument("--out-saturation", required=True, help="Output saturation CSV")
    parser.add_argument("--out-plot", default=None, help="Output saturation plot PNG")
    return parser.parse_args()


def alignment_stats(records):
    """Compute basic alignment statistics."""
    if not records:
        return {"n_seqs": 0, "aln_length": 0, "gap_pct": 0, "informative_sites": 0}

    seqs = [str(r.seq).upper() for r in records]
    n_seqs = len(seqs)
    aln_length = len(seqs[0]) if seqs else 0

    # Gap percentage
    total_chars = sum(len(s) for s in seqs)
    total_gaps = sum(s.count("-") for s in seqs)
    gap_pct = (total_gaps / total_chars * 100) if total_chars > 0 else 0

    # Parsimony-informative sites: sites with ≥2 different chars each appearing ≥2 times
    informative = 0
    for col in range(aln_length):
        chars = [s[col] for s in seqs if s[col] not in ("-", "N")]
        counts = Counter(chars)
        n_informative_states = sum(1 for c in counts.values() if c >= 2)
        if n_informative_states >= 2:
            informative += 1

    return {
        "n_seqs": n_seqs,
        "aln_length": aln_length,
        "gap_pct": round(gap_pct, 2),
        "informative_sites": informative,
    }


def count_ts_tv_by_position(seqs, n_pairs_max=500):
    """
    Count transitions and transversions per codon position
    across pairwise sequence comparisons.
    Returns dict: {codon_pos: {"ts": [...], "tv": [...]}}
    """
    results = {1: {"ts": [], "tv": []}, 2: {"ts": [], "tv": []}, 3: {"ts": [], "tv": []}}
    aln_len = len(seqs[0])

    # Subsample pairs if too many sequences
    pairs = list(combinations(range(len(seqs)), 2))
    if len(pairs) > n_pairs_max:
        rng = np.random.default_rng(42)
        indices = rng.choice(len(pairs), size=n_pairs_max, replace=False)
        pairs = [pairs[i] for i in indices]

    for i, j in pairs:
        ts_counts = {1: 0, 2: 0, 3: 0}
        tv_counts = {1: 0, 2: 0, 3: 0}

        for pos in range(aln_len):
            a, b = seqs[i][pos], seqs[j][pos]
            if a in ("-", "N") or b in ("-", "N"):
                continue
            if a == b:
                continue
            codon_pos = (pos % 3) + 1  # 1, 2, or 3
            if (a, b) in TRANSITIONS:
                ts_counts[codon_pos] += 1
            else:
                tv_counts[codon_pos] += 1

        for cp in (1, 2, 3):
            results[cp]["ts"].append(ts_counts[cp])
            results[cp]["tv"].append(tv_counts[cp])

    return results


def xia_saturation_index(ts_tv_data):
    """
    Simplified Xia's Iss (Index of Substitution Saturation).

    Concept: If Ts/Tv ratio approaches 0.5 (random), sequences are saturated.
    Iss = observed entropy / critical entropy.
    If Iss >= Iss_c → saturated.

    Simplified version: Iss ≈ proportion of Ts among all substitutions.
    Iss_c ≈ 0.7 (empirical threshold for coding sequences).
    """
    results = {}
    for cp in (1, 2, 3):
        ts_total = sum(ts_tv_data[cp]["ts"])
        tv_total = sum(ts_tv_data[cp]["tv"])
        total = ts_total + tv_total

        if total == 0:
            iss = 0.0
        else:
            # Iss approximation: when Ts/(Ts+Tv) drops toward 0.5, saturation
            ts_ratio = ts_total / total
            # Invert: if ts_ratio is around 0.5, high saturation
            iss = 1.0 - abs(ts_ratio - 0.5) * 2  # 0 = no saturation, 1 = full

        iss_c = 0.7  # Critical threshold
        saturated = iss >= iss_c

        results[cp] = {
            "iss": round(iss, 4),
            "iss_c": iss_c,
            "ts_total": ts_total,
            "tv_total": tv_total,
            "saturated": saturated,
        }

    return results


def plot_saturation(ts_tv_data, out_path, gene_name):
    """Plot Ts vs Tv scatter for each codon position."""
    if not HAS_MPL:
        print("[warn] matplotlib not available, skipping plot")
        return

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.suptitle(f"{gene_name} — Substitution Saturation", fontsize=13)

    for idx, cp in enumerate([1, 2, 3]):
        ax = axes[idx]
        ts = ts_tv_data[cp]["ts"]
        tv = ts_tv_data[cp]["tv"]

        ax.scatter(tv, ts, alpha=0.4, s=15, color=["#2E86AB", "#E8430C", "#44AF69"][idx])
        max_val = max(max(ts, default=1), max(tv, default=1))
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="Ts=Tv (saturation)")
        ax.set_xlabel("Transversions")
        ax.set_ylabel("Transitions")
        ax.set_title(f"Codon position {cp}")
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[plot] Saved: {out_path}")


def main():
    args = parse_args()

    # Read alignment
    records = list(SeqIO.parse(args.input, "fasta"))
    seqs = [str(r.seq).upper() for r in records]

    # Alignment stats
    stats = alignment_stats(records)
    with open(args.out_stats, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "n_seqs", "aln_length", "gap_pct", "informative_sites"])
        w.writerow([args.gene, stats["n_seqs"], stats["aln_length"],
                     stats["gap_pct"], stats["informative_sites"]])
    print(f"[{args.gene}] Alignment: {stats['n_seqs']} seqs, {stats['aln_length']} cols, "
          f"{stats['gap_pct']}% gaps, {stats['informative_sites']} informative sites")

    # Saturation test
    if len(seqs) < 4:
        print(f"[{args.gene}] Too few sequences for saturation test")
        with open(args.out_saturation, "w") as f:
            f.write("gene,codon_pos,iss,iss_c,ts_total,tv_total,saturated\n")
        return

    ts_tv = count_ts_tv_by_position(seqs)
    sat_results = xia_saturation_index(ts_tv)

    with open(args.out_saturation, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "codon_pos", "iss", "iss_c", "ts_total", "tv_total", "saturated"])
        for cp in (1, 2, 3):
            r = sat_results[cp]
            w.writerow([args.gene, cp, r["iss"], r["iss_c"],
                        r["ts_total"], r["tv_total"], r["saturated"]])
            status = "⚠️ SATURATED" if r["saturated"] else "✅ OK"
            print(f"[{args.gene}] Codon pos {cp}: Iss={r['iss']:.4f} "
                  f"(Iss_c={r['iss_c']}) Ts={r['ts_total']} Tv={r['tv_total']} → {status}")

    # Plot
    if args.out_plot:
        plot_saturation(ts_tv, args.out_plot, args.gene)


if __name__ == "__main__":
    main()
