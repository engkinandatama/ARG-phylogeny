#!/usr/bin/env python3
"""
Build strain-gene co-occurrence network and co-selection analysis.

Usage:
  build_network.py --input-dir . --out-gexf network.gexf \
                   --out-coselection matrix.csv --out-heatmap heatmap.png
"""

import argparse
import csv
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd
import networkx as nx

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False


def parse_args():
    parser = argparse.ArgumentParser(description="Build ARG network")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--out-gexf", required=True)
    parser.add_argument("--out-coselection", required=True)
    parser.add_argument("--out-heatmap", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)

    # Build network
    G = nx.Graph()
    gene_strains = defaultdict(set)

    for meta in sorted(input_dir.glob("*_meta.csv")):
        gene = meta.stem.replace("_meta", "")
        try:
            df = pd.read_csv(meta)
        except Exception:
            continue

        for _, row in df.iterrows():
            strain = str(row.get("acc", ""))
            country = str(row.get("country", "NA")).split(":")[0] if pd.notna(row.get("country")) else "NA"
            date = str(row.get("date", "NA"))

            G.add_node(strain, type="strain", country=country, date=date)
            G.add_node(gene, type="gene")
            G.add_edge(strain, gene)
            gene_strains[gene].add(strain)

    nx.write_gexf(G, args.out_gexf)
    print(f"[network] Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Co-selection analysis (Fisher's exact test)
    all_strains = set()
    for s in gene_strains.values():
        all_strains |= s

    co_results = []
    genes = sorted(gene_strains.keys())
    for g1, g2 in combinations(genes, 2):
        s1 = gene_strains[g1]
        s2 = gene_strains[g2]
        both = len(s1 & s2)
        only1 = len(s1 - s2)
        only2 = len(s2 - s1)
        neither = len(all_strains - s1 - s2)

        jaccard = both / len(s1 | s2) if (s1 | s2) else 0

        pval = 1.0
        if HAS_SCIPY:
            table = [[both, only1], [only2, neither]]
            _, pval = scipy_stats.fisher_exact(table)

        co_results.append({
            "gene_a": g1, "gene_b": g2,
            "co_occur": both, "jaccard": round(jaccard, 4),
            "p_value": round(pval, 6)
        })

    with open(args.out_coselection, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_a", "gene_b", "co_occur", "jaccard", "p_value"])
        w.writeheader()
        w.writerows(co_results)

    # Heatmap
    if args.out_heatmap and HAS_PLOT and genes:
        matrix = pd.DataFrame(0.0, index=genes, columns=genes)
        for r in co_results:
            matrix.loc[r["gene_a"], r["gene_b"]] = r["jaccard"]
            matrix.loc[r["gene_b"], r["gene_a"]] = r["jaccard"]
        for g in genes:
            matrix.loc[g, g] = 1.0

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(matrix, annot=True, cmap="YlOrRd", vmin=0, vmax=1,
                    fmt=".2f", ax=ax, square=True)
        ax.set_title("ARG Co-occurrence (Jaccard Index)")
        plt.tight_layout()
        plt.savefig(args.out_heatmap, dpi=150)
        plt.close()
        print(f"[network] Heatmap: {args.out_heatmap}")


if __name__ == "__main__":
    main()
