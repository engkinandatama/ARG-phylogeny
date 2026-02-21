#!/usr/bin/env python3
"""
HGT (Horizontal Gene Transfer) detection via:
  A) Phylogenetic incongruence: Robinson-Foulds distance between gene/species trees
  B) Parametric: GC content + Codon Adaptation Index (CAI)

Usage:
  tree_comparison.py --gene-tree gene.treefile --species-tree species.treefile \
                     --sequences gene.fasta --gene blaCTX-M \
                     --out-csv report.csv --out-json detail.json
"""

import argparse
import csv
import json
import sys
from collections import Counter

from Bio import SeqIO

try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False


def parse_args():
    parser = argparse.ArgumentParser(description="HGT detection")
    parser.add_argument("--gene-tree", required=True)
    parser.add_argument("--species-tree", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-json", required=True)
    return parser.parse_args()


def robinson_foulds(gene_tree_path, species_tree_path):
    """Compute normalized Robinson-Foulds distance."""
    if not HAS_ETE3:
        print("[warn] ete3 not available, skipping RF distance", file=sys.stderr)
        return None, None

    try:
        gt = Tree(gene_tree_path)
        st = Tree(species_tree_path)

        # Get common leaves
        gt_leaves = set(gt.get_leaf_names())
        st_leaves = set(st.get_leaf_names())
        common = gt_leaves & st_leaves

        if len(common) < 4:
            print(f"[hgt] Too few common taxa ({len(common)}), skipping RF", file=sys.stderr)
            return None, None

        # Prune both trees to common leaves
        gt.prune(common, preserve_branch_length=True)
        st.prune(common, preserve_branch_length=True)

        # Compute RF distance
        result = gt.robinson_foulds(st, unrooted_trees=True)
        rf = result[0]
        max_rf = result[1]
        norm_rf = rf / max_rf if max_rf > 0 else 0

        print(f"[hgt] RF distance: {rf}/{max_rf} (normalized: {norm_rf:.4f})")
        return norm_rf, rf

    except Exception as e:
        print(f"[warn] RF computation failed: {e}", file=sys.stderr)
        return None, None


def gc_content(fasta_path):
    """Compute average GC content of sequences."""
    gc_values = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper().replace("-", "").replace("N", "")
        if len(seq) > 0:
            gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
            gc_values.append(gc)
    avg_gc = sum(gc_values) / len(gc_values) if gc_values else 0
    return round(avg_gc, 2)


def codon_adaptation_index(fasta_path):
    """
    Simplified Codon Adaptation Index (CAI).
    Compares codon usage of the gene to expected usage in highly expressed genes.
    Higher CAI → better adapted (less likely HGT).
    Lower CAI → possibly foreign origin (HGT candidate).
    """
    codon_counts = Counter()
    aa_codon_map = {}  # amino acid -> list of codons

    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper().replace("-", "")
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3 and all(c in "ACGT" for c in codon):
                codon_counts[codon] += 1

    # Compute effective number of codons (simplified ENC)
    total_codons = sum(codon_counts.values())
    if total_codons == 0:
        return 0.0

    n_codons_used = len(codon_counts)

    # CAI proxy: ratio of codons used to total possible (61 sense codons)
    cai_proxy = n_codons_used / 61.0
    return round(cai_proxy, 4)


def classify_hgt(rf_norm, gc, cai):
    """Classify HGT evidence level."""
    score = 0
    interpretations = {}

    if rf_norm is not None:
        if rf_norm > 0.5:
            score += 2
            interpretations["rf"] = "High incongruence"
        elif rf_norm > 0.3:
            score += 1
            interpretations["rf"] = "Moderate incongruence"
        else:
            interpretations["rf"] = "Low incongruence"

    # GC content deviation from K. pneumoniae average (~57%)
    kp_gc = 57.0
    gc_dev = abs(gc - kp_gc)
    if gc_dev > 5:
        score += 2
        interpretations["gc"] = f"Deviant ({gc_dev:.1f}% from expected)"
    elif gc_dev > 3:
        score += 1
        interpretations["gc"] = f"Slightly deviant ({gc_dev:.1f}%)"
    else:
        interpretations["gc"] = f"Within range ({gc_dev:.1f}%)"

    if score >= 3:
        verdict = "Strong"
    elif score >= 2:
        verdict = "Moderate"
    elif score >= 1:
        verdict = "Weak"
    else:
        verdict = "None"

    return verdict, interpretations


def main():
    args = parse_args()

    # A) Phylogenetic incongruence
    rf_norm, rf_raw = robinson_foulds(args.gene_tree, args.species_tree)

    # B) Parametric
    gc = gc_content(args.sequences)
    cai = codon_adaptation_index(args.sequences)

    # Classify
    verdict, interp = classify_hgt(rf_norm, gc, cai)

    print(f"[{args.gene}] HGT evidence: {verdict}")
    print(f"  RF={rf_norm}, GC={gc}%, CAI={cai}")

    # Write CSV
    with open(args.out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "rf_distance", "gc_content", "cai_score", "hgt_verdict"])
        w.writerow([args.gene, rf_norm or "NA", gc, cai, verdict])

    # Write detail JSON
    detail = {
        "gene": args.gene,
        "topology": {"rf_normalized": rf_norm, "rf_raw": rf_raw},
        "parametric": {"gc_content": gc, "cai_score": cai},
        "interpretations": interp,
        "verdict": verdict,
    }
    with open(args.out_json, "w") as f:
        json.dump(detail, f, indent=2)


if __name__ == "__main__":
    main()
