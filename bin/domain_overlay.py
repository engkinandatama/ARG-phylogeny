#!/usr/bin/env python3
"""
Domain annotation (via HMMER/Pfam) and domain-selection overlay plotting.

Two modes:
  annotate: Run hmmscan and parse domain hits to JSON
  plot:     Overlay selected sites on domain architecture (PNG + HTML)

Usage:
  domain_overlay.py annotate --input gene.fasta --gene blaCTX-M --output domains.json
  domain_overlay.py plot --domains domains.json --selection sig.csv \
                         --out-png overlay.png --out-html overlay.html
"""

import argparse
import csv
import json
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def parse_args():
    parser = argparse.ArgumentParser(description="Domain annotation & overlay")
    subparsers = parser.add_subparsers(dest="command")

    # Annotate
    ann = subparsers.add_parser("annotate")
    ann.add_argument("--input", required=True)
    ann.add_argument("--gene", required=True)
    ann.add_argument("--output", required=True)
    ann.add_argument("--pfam-db", default="Pfam-A.hmm")

    # Plot
    plt_parser = subparsers.add_parser("plot")
    plt_parser.add_argument("--domains", required=True)
    plt_parser.add_argument("--selection", required=True)
    plt_parser.add_argument("--gene", default="gene")
    plt_parser.add_argument("--out-png", required=True)
    plt_parser.add_argument("--out-html", required=True)

    return parser.parse_args()


def annotate_domains(input_fasta, gene_name, output_json, pfam_db="Pfam-A.hmm"):
    """Run hmmscan and parse domain hits."""
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        json.dump({"gene": gene_name, "domains": []}, open(output_json, "w"), indent=2)
        return

    # Take first/representative sequence, translate
    rep = records[0]
    seq = str(rep.seq).upper().replace("-", "")
    if len(seq) % 3 != 0:
        seq = seq[:len(seq) - len(seq) % 3]
    aa_seq = str(Seq(seq).translate(table=11)).rstrip("*")

    # Write AA to temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
        f.write(f">{rep.id}\n{aa_seq}\n")
        aa_tmp = f.name

    domains = []
    try:
        # Run hmmscan
        with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as tbl:
            tbl_path = tbl.name

        cmd = ["hmmscan", "--domtblout", tbl_path, "--noali", "-E", "1e-5",
               pfam_db, aa_tmp]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            # Parse domtblout
            with open(tbl_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 23:
                        domains.append({
                            "name": parts[0],
                            "accession": parts[1],
                            "description": " ".join(parts[22:]),
                            "start": int(parts[17]),  # ali_from
                            "end": int(parts[18]),     # ali_to
                            "evalue": float(parts[6]),
                        })
            print(f"[domain] {gene_name}: {len(domains)} domains found")
        else:
            print(f"[domain] hmmscan not available or failed, using empty domains")

    except FileNotFoundError:
        print(f"[domain] hmmscan not found, skipping domain annotation")

    result = {"gene": gene_name, "seq_length_aa": len(aa_seq), "domains": domains}
    with open(output_json, "w") as f:
        json.dump(result, f, indent=2)


def plot_overlay(domains_json, selection_csv, gene_name, out_png, out_html):
    """Plot domain architecture with selection site overlay."""
    with open(domains_json) as f:
        dom_data = json.load(f)

    domains = dom_data.get("domains", [])
    seq_len = dom_data.get("seq_length_aa", 300)

    # Parse selection sites
    pos_sites = []
    neg_sites = []
    try:
        with open(selection_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                site = int(row.get("site", 0))
                if site > 0:
                    pos_sites.append(site)
    except Exception:
        pass

    # Static plot (matplotlib)
    if HAS_MPL:
        fig, ax = plt.subplots(figsize=(12, 3))
        colors = plt.cm.Set3.colors

        # Draw protein backbone
        ax.barh(0, seq_len, height=0.3, color="#E0E0E0", edgecolor="none")

        # Draw domains
        for i, d in enumerate(domains):
            color = colors[i % len(colors)]
            width = d["end"] - d["start"]
            rect = mpatches.FancyBboxPatch(
                (d["start"], -0.2), width, 0.4,
                boxstyle="round,pad=0.02", facecolor=color,
                edgecolor="black", linewidth=0.5
            )
            ax.add_patch(rect)
            ax.text(d["start"] + width/2, 0, d["name"], ha="center", va="center", fontsize=7)

        # Mark selected sites
        for site in pos_sites:
            if site <= seq_len:
                ax.plot(site, 0.35, "v", color="red", markersize=6, zorder=5)

        ax.set_xlim(0, seq_len)
        ax.set_ylim(-0.5, 0.6)
        ax.set_xlabel("Amino acid position")
        ax.set_title(f"{gene_name} — Domain Architecture + Selection Sites")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()
        print(f"[plot] Static: {out_png}")

    # Interactive plot (Plotly)
    if HAS_PLOTLY:
        fig = go.Figure()

        # Backbone
        fig.add_shape(type="rect", x0=0, x1=seq_len, y0=-0.15, y1=0.15,
                       fillcolor="#E0E0E0", line_width=0)

        # Domains
        colors_plotly = ["#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                          "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"]
        for i, d in enumerate(domains):
            fig.add_shape(type="rect",
                          x0=d["start"], x1=d["end"], y0=-0.2, y1=0.2,
                          fillcolor=colors_plotly[i % len(colors_plotly)],
                          line=dict(width=1))
            fig.add_annotation(x=(d["start"]+d["end"])/2, y=0,
                              text=d["name"], showarrow=False, font_size=10)

        # Selected sites
        if pos_sites:
            fig.add_trace(go.Scatter(
                x=pos_sites, y=[0.3]*len(pos_sites),
                mode="markers", marker=dict(symbol="triangle-down", size=10, color="red"),
                name="Positively selected sites",
                hovertemplate="Site %{x}<extra>Positive selection</extra>"
            ))

        fig.update_layout(
            title=f"{gene_name} — Domain Architecture + Selection Sites",
            xaxis_title="Amino acid position",
            yaxis=dict(visible=False, range=[-0.5, 0.6]),
            showlegend=True, height=250
        )
        fig.write_html(out_html)
        print(f"[plot] Interactive: {out_html}")


def main():
    args = parse_args()
    if args.command == "annotate":
        annotate_domains(args.input, args.gene, args.output,
                          pfam_db=getattr(args, "pfam_db", "Pfam-A.hmm"))
    elif args.command == "plot":
        plot_overlay(args.domains, args.selection, args.gene, args.out_png, args.out_html)
    else:
        print("Usage: domain_overlay.py {annotate|plot}")
        sys.exit(1)


if __name__ == "__main__":
    main()
