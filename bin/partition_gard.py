#!/usr/bin/env python3
"""
Parse GARD JSON results and split alignment into partitions.

If GARD detected breakpoints, the alignment is split at those positions.
If no breakpoints, the full alignment is copied as partition_1.

Usage:
  partition_gard.py --gard-json gene_GARD.json --alignment gene.trimmed.fasta \
                    --gene blaCTX-M --outdir partitions/
"""

import argparse
import json
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description="Parse GARD and partition alignment")
    parser.add_argument("--gard-json", required=True, help="GARD JSON output")
    parser.add_argument("--alignment", required=True, help="Input alignment FASTA")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--outdir", required=True, help="Output directory for partitions")
    return parser.parse_args()


def extract_breakpoints(gard_json_path):
    """Extract breakpoint positions from GARD JSON."""
    try:
        with open(gard_json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"[warn] Could not parse GARD JSON: {e}", file=sys.stderr)
        return []

    breakpoints = []

    # GARD JSON structure varies by HyPhy version
    # Try common locations
    if "breakpoints" in data:
        bp = data["breakpoints"]
        if isinstance(bp, list):
            breakpoints = sorted(bp)
        elif isinstance(bp, dict):
            # Some versions use nested structure
            for key, val in bp.items():
                if isinstance(val, (int, float)):
                    breakpoints.append(int(val))
    elif "improvements" in data:
        # Alternative structure
        for imp in data.get("improvements", {}).values():
            if "breakpoint" in imp:
                breakpoints.append(int(imp["breakpoint"]))
        breakpoints = sorted(breakpoints)

    return breakpoints


def split_alignment(records, breakpoints, gene_name, outdir):
    """Split alignment at breakpoint positions."""
    if not records:
        return

    aln_len = len(str(records[0].seq))

    if not breakpoints:
        # No breakpoints — write full alignment as single partition
        out_path = f"{outdir}/{gene_name}_partition_1.fasta"
        SeqIO.write(records, out_path, "fasta")
        print(f"[partition] No breakpoints detected. Full alignment as partition 1.")
        return

    # Add boundaries
    boundaries = [0] + breakpoints + [aln_len]
    print(f"[partition] Breakpoints: {breakpoints}")
    print(f"[partition] Splitting into {len(boundaries)-1} partitions")

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]

        partition_records = []
        for rec in records:
            seq_slice = str(rec.seq)[start:end]
            partition_records.append(
                SeqRecord(Seq(seq_slice), id=rec.id, description="")
            )

        out_path = f"{outdir}/{gene_name}_partition_{i+1}.fasta"
        SeqIO.write(partition_records, out_path, "fasta")
        print(f"[partition] Partition {i+1}: positions {start+1}–{end} ({end-start} nt)")


def main():
    args = parse_args()

    # Parse breakpoints
    breakpoints = extract_breakpoints(args.gard_json)

    # Read alignment
    records = list(SeqIO.parse(args.alignment, "fasta"))

    # Split
    split_alignment(records, breakpoints, args.gene, args.outdir)


if __name__ == "__main__":
    main()
