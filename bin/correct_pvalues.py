#!/usr/bin/env python3
"""
Correct p-values from HyPhy selection analysis for multiple testing.

Applies Benjamini-Hochberg (FDR) or Holm-Bonferroni correction across
all sites/methods for a given gene.

Usage:
  correct_pvalues.py --input-dir . --gene blaCTX-M --method BH \
                     --out-summary summary.csv --out-significant sig.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    # Fallback: implement simple BH
    multipletests = None


def parse_args():
    parser = argparse.ArgumentParser(description="Multiple testing correction for HyPhy")
    parser.add_argument("--input-dir", required=True, help="Dir with HyPhy JSON files")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--method", default="BH", choices=["BH", "holm", "bonferroni"])
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-significant", required=True)
    return parser.parse_args()


def parse_hyphy_pvalues(json_path, method_name):
    """Extract site-level p-values from HyPhy JSON."""
    try:
        with open(json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        return []

    sites = []

    # Site-level methods: SLAC, FEL, MEME
    if "MLE" in data:
        mle = data["MLE"]
        headers = mle.get("headers", [])
        content = mle.get("content", {})

        # Find p-value column
        p_col = None
        for i, h in enumerate(headers):
            if isinstance(h, list) and len(h) > 0:
                label = h[0].lower()
            elif isinstance(h, str):
                label = h.lower()
            else:
                continue
            if "p-value" in label or label == "p":
                p_col = i
                break

        if p_col is not None:
            # content can be dict (site_idx -> row) or list
            if isinstance(content, dict):
                for site_key, row in content.items():
                    if isinstance(row, list) and len(row) > p_col:
                        try:
                            pval = float(row[p_col])
                            sites.append({"site": site_key, "p_value": pval})
                        except (ValueError, TypeError):
                            pass
            elif isinstance(content, list):
                for i, row in enumerate(content):
                    if isinstance(row, list) and len(row) > p_col:
                        try:
                            pval = float(row[p_col])
                            sites.append({"site": str(i+1), "p_value": pval})
                        except (ValueError, TypeError):
                            pass

    # BUSTED: gene-wide p-value
    elif "test results" in data:
        pval = data["test results"].get("p-value")
        if pval is not None:
            sites.append({"site": "gene-wide", "p_value": float(pval)})

    return sites


def apply_correction(all_pvalues, method="BH"):
    """Apply multiple testing correction."""
    if not all_pvalues:
        return []

    pvals = np.array([x["p_value"] for x in all_pvalues])

    if multipletests is not None:
        method_map = {"BH": "fdr_bh", "holm": "holm", "bonferroni": "bonferroni"}
        reject, corrected, _, _ = multipletests(pvals, method=method_map.get(method, "fdr_bh"))
        for i, entry in enumerate(all_pvalues):
            entry["q_value"] = float(corrected[i])
            entry["reject"] = bool(reject[i])
    else:
        # Simple BH fallback
        n = len(pvals)
        sorted_idx = np.argsort(pvals)
        for rank, idx in enumerate(sorted_idx, 1):
            q = pvals[idx] * n / rank
            all_pvalues[idx]["q_value"] = min(q, 1.0)
            all_pvalues[idx]["reject"] = all_pvalues[idx]["q_value"] < 0.05

    return all_pvalues


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)

    # Collect p-values from all methods
    all_entries = []
    methods_found = []

    for json_file in sorted(input_dir.glob(f"{args.gene}_*.json")):
        method_name = json_file.stem.replace(f"{args.gene}_", "")
        sites = parse_hyphy_pvalues(json_file, method_name)
        for s in sites:
            s["method"] = method_name
            s["gene"] = args.gene
        all_entries.extend(sites)
        if sites:
            methods_found.append(method_name)

    print(f"[correction] {args.gene}: {len(all_entries)} p-values from {methods_found}")

    # Apply correction
    corrected = apply_correction(all_entries, method=args.method)

    # Write summary
    with open(args.out_summary, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "method", "site", "p_value", "q_value", "significant"])
        for entry in corrected:
            w.writerow([
                entry["gene"], entry["method"], entry["site"],
                f"{entry['p_value']:.6f}", f"{entry['q_value']:.6f}",
                entry.get("reject", False)
            ])

    # Write significant sites only
    sig_entries = [e for e in corrected if e.get("reject", False)]
    with open(args.out_significant, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "method", "site", "p_value", "q_value"])
        for entry in sig_entries:
            w.writerow([
                entry["gene"], entry["method"], entry["site"],
                f"{entry['p_value']:.6f}", f"{entry['q_value']:.6f}"
            ])

    print(f"[correction] {len(sig_entries)} significant after {args.method} correction "
          f"(threshold={args.p_threshold})")


if __name__ == "__main__":
    main()
