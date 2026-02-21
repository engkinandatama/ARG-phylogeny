#!/usr/bin/env python3
"""
Generate per-gene markdown report from analysis results using Jinja2 template.

Usage:
  generate_report.py --gene blaCTX-M --qc-stats qc.csv --aln-stats aln.csv \
                     --saturation sat.csv --selection sel.csv --hgt hgt.csv \
                     --domains domains.json --temporal temporal.csv \
                     --template report.jinja2 --out-md report.md
"""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

from jinja2 import Template


def parse_args():
    parser = argparse.ArgumentParser(description="Generate analysis report")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--qc-stats", required=True)
    parser.add_argument("--aln-stats", required=True)
    parser.add_argument("--saturation", required=True)
    parser.add_argument("--selection", required=True)
    parser.add_argument("--hgt", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument("--temporal", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-html", default=None)
    return parser.parse_args()


def read_csv_first_row(path):
    """Read first data row of a CSV as dict."""
    try:
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                return row
    except Exception:
        return {}
    return {}


def main():
    args = parse_args()

    # Load data
    qc = read_csv_first_row(args.qc_stats)
    aln = read_csv_first_row(args.aln_stats)
    sat = read_csv_first_row(args.saturation)
    sel = read_csv_first_row(args.selection)
    hgt = read_csv_first_row(args.hgt)
    temporal = read_csv_first_row(args.temporal)

    try:
        with open(args.domains) as f:
            domains_data = json.load(f)
    except Exception:
        domains_data = {"domains": []}

    # Load template
    template_text = Path(args.template).read_text()
    template = Template(template_text)

    # Render
    context = {
        "gene": args.gene,
        "organism": "Klebsiella pneumoniae",
        "version": "0.1.0",
        "date": datetime.now().strftime("%Y-%m-%d"),
        "n_raw": qc.get("input", "NA"),
        "n_clustered": "NA",
        "n_clean": qc.get("kept", "NA"),
        "aln_length": aln.get("aln_length", "NA"),
        "aln_trimmed_length": aln.get("aln_length", "NA"),
        "gap_pct": aln.get("gap_pct", "NA"),
        "informative_sites": aln.get("informative_sites", "NA"),
        "iss_value": sat.get("iss", "NA"),
        "iss_critical": sat.get("iss_c", "NA"),
        "saturation_status": "SATURATED" if sat.get("saturated") == "True" else "OK",
        "tree_model": "MFP (auto)",
        "tree_lnl": "see log",
        "pct_supported": "NA",
        "rooting_method": "midpoint",
        "gard_breakpoints": 0,
        "gard_positions": "",
        "correction_method": "Benjamini-Hochberg",
        "p_thresh": 0.05,
        "q_thresh": 0.1,
        "busted_p": sel.get("p_value", "NA") if sel.get("method") == "BUSTED" else "NA",
        "slac_pos": "NA", "slac_neg": "NA",
        "fel_pos": "NA", "fel_neg": "NA",
        "meme_pos": "NA",
        "fubar_pos": "NA", "fubar_neg": "NA",
        "absrel_branches": "NA",
        "relax_k": "NA", "relax_p": "NA",
        "rf_distance": hgt.get("rf_distance", "NA"),
        "rf_interpretation": "",
        "au_pvalue": "NA",
        "au_interpretation": "",
        "gc_gene": hgt.get("gc_content", "NA"),
        "gc_interpretation": "",
        "cai_score": hgt.get("cai_score", "NA"),
        "cai_interpretation": "",
        "hgt_verdict": hgt.get("hgt_verdict", "NA"),
        "domains": domains_data.get("domains", []),
        "variant_names": [],
        "temporal_signal": temporal.get("r_squared") not in (None, "NA", ""),
        "temporal_r2": temporal.get("r_squared", "NA"),
        "temporal_slope": temporal.get("slope", "NA"),
        "temporal_date_range": temporal.get("date_range", "NA"),
    }

    report_md = template.render(**context)

    with open(args.out_md, "w") as f:
        f.write(report_md)
    print(f"[report] Generated: {args.out_md}")

    if args.out_html:
        try:
            import markdown
            html = markdown.markdown(report_md, extensions=["tables"])
            with open(args.out_html, "w") as f:
                f.write(f"<html><head><style>body{{font-family:sans-serif;max-width:800px;margin:auto;padding:20px}}table{{border-collapse:collapse}}td,th{{border:1px solid #ddd;padding:6px}}</style></head><body>{html}</body></html>")
            print(f"[report] HTML: {args.out_html}")
        except ImportError:
            print("[report] python-markdown not available, skipping HTML")


if __name__ == "__main__":
    main()
