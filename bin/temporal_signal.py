#!/usr/bin/env python3
"""
Temporal signal analysis: root-to-tip regression.
Correlates root-to-tip distance with collection date to detect clock-like signal.

Usage:
  temporal_signal.py --tree gene.treefile --meta gene_meta.csv \
                     --gene blaCTX-M --out-csv temporal.csv --out-plot temporal.png
"""

import argparse
import csv
import re
import sys

import numpy as np

try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def parse_args():
    parser = argparse.ArgumentParser(description="Temporal signal analysis")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-plot", default=None)
    return parser.parse_args()


def extract_year(date_str):
    """Extract year from various date formats."""
    if not date_str or date_str == "NA":
        return None
    # Try patterns: YYYY, DD-MMM-YYYY, YYYY-MM-DD
    match = re.search(r"(\d{4})", str(date_str))
    if match:
        year = int(match.group(1))
        if 1980 <= year <= 2030:
            return year
    return None


def main():
    args = parse_args()

    if not HAS_ETE3:
        print("[temporal] ete3 not available, skipping")
        with open(args.out_csv, "w") as f:
            f.write("gene,r_squared,slope,intercept,n_samples,date_range\n")
            f.write(f"{args.gene},NA,NA,NA,0,NA\n")
        return

    # Load tree
    try:
        t = Tree(args.tree)
    except Exception as e:
        print(f"[temporal] Cannot load tree: {e}", file=sys.stderr)
        with open(args.out_csv, "w") as f:
            f.write("gene,r_squared,slope,intercept,n_samples,date_range\n")
            f.write(f"{args.gene},NA,NA,NA,0,NA\n")
        return

    # Load metadata dates
    dates = {}
    try:
        with open(args.meta) as f:
            reader = csv.DictReader(f)
            for row in reader:
                acc = row.get("acc", "")
                date = row.get("date", "")
                year = extract_year(date)
                if year:
                    dates[acc] = year
    except Exception:
        pass

    # Get root-to-tip distances
    root = t.get_tree_root()
    data_points = []
    for leaf in t.get_leaves():
        if leaf.name in dates:
            dist = t.get_distance(root, leaf)
            data_points.append((dates[leaf.name], dist))

    if len(data_points) < 5:
        print(f"[temporal] Too few dated samples ({len(data_points)})")
        with open(args.out_csv, "w") as f:
            f.write("gene,r_squared,slope,intercept,n_samples,date_range\n")
            f.write(f"{args.gene},NA,NA,NA,{len(data_points)},NA\n")
        return

    years = np.array([p[0] for p in data_points])
    dists = np.array([p[1] for p in data_points])

    # Linear regression
    slope, intercept = np.polyfit(years, dists, 1)
    ss_res = np.sum((dists - (slope * years + intercept)) ** 2)
    ss_tot = np.sum((dists - np.mean(dists)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    date_range = f"{int(years.min())}-{int(years.max())}"
    print(f"[{args.gene}] Temporal: R²={r_squared:.4f}, slope={slope:.6f}, range={date_range}")

    # Write CSV
    with open(args.out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "r_squared", "slope", "intercept", "n_samples", "date_range"])
        w.writerow([args.gene, f"{r_squared:.4f}", f"{slope:.6f}",
                     f"{intercept:.6f}", len(data_points), date_range])

    # Plot
    if args.out_plot and HAS_MPL:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(years, dists, alpha=0.5, s=20, color="#2E86AB")
        x_line = np.linspace(years.min(), years.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, "r--", linewidth=1.5,
                label=f"R²={r_squared:.4f}")
        ax.set_xlabel("Collection Year")
        ax.set_ylabel("Root-to-tip Distance")
        ax.set_title(f"{args.gene} — Temporal Signal")
        ax.legend()
        plt.tight_layout()
        plt.savefig(args.out_plot, dpi=150)
        plt.close()
        print(f"[plot] Saved: {args.out_plot}")


if __name__ == "__main__":
    main()
