#!/usr/bin/env python3
"""
QC / Filter / Dedup sequences for ARG phylogenomics.

Filters:
  1. Minimum nucleotide length
  2. Maximum ambiguity percentage
  3. Reading frame check (length divisible by 3)
  4. Internal stop codon removal
  5. Exact sequence deduplication

Outputs:
  - Cleaned FASTA file
  - QC statistics CSV

Usage:
  qc_sequences.py --input gene.fasta --output gene.clean.fasta \
                   --stats gene_qc_stats.csv --gene blaCTX-M
"""

import argparse
import csv
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Standard stop codons (universal + bacterial)
STOP_CODONS = {"TAA", "TAG", "TGA"}
AMBIGUOUS_BASES = set("NRYWSKMBDHV")


def parse_args():
    parser = argparse.ArgumentParser(description="QC filter sequences")
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output clean FASTA")
    parser.add_argument("--stats", required=True, help="Output QC stats CSV")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--min-len", type=int, default=300, help="Min length (nt)")
    parser.add_argument("--max-ambig", type=float, default=5.0, help="Max ambiguity %%")
    return parser.parse_args()


def qc_filter(input_fasta, min_len, max_ambig):
    """Filter sequences and return kept sequences + statistics."""
    stats = {
        "input": 0,
        "removed_short": 0,
        "removed_ambig": 0,
        "removed_frame": 0,
        "removed_stop": 0,
        "removed_dup": 0,
        "kept": 0,
    }

    kept = []
    seen_seqs = set()

    for rec in SeqIO.parse(input_fasta, "fasta"):
        stats["input"] += 1
        seq = str(rec.seq).upper().replace("U", "T")

        # Filter 1: Minimum length
        if len(seq) < min_len:
            stats["removed_short"] += 1
            continue

        # Filter 2: Maximum ambiguity
        n_ambig = sum(1 for c in seq if c in AMBIGUOUS_BASES)
        ambig_pct = (n_ambig / len(seq)) * 100
        if ambig_pct > max_ambig:
            stats["removed_ambig"] += 1
            continue

        # Filter 3: Reading frame (divisible by 3)
        if len(seq) % 3 != 0:
            stats["removed_frame"] += 1
            continue

        # Filter 4: Internal stop codons
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        # Check all codons except the last one (which can be a stop)
        if any(c in STOP_CODONS for c in codons[:-1]):
            stats["removed_stop"] += 1
            continue

        # Filter 5: Exact sequence deduplication
        if seq in seen_seqs:
            stats["removed_dup"] += 1
            continue

        seen_seqs.add(seq)
        kept.append(
            SeqRecord(Seq(seq), id=rec.id, description=rec.description)
        )

    stats["kept"] = len(kept)
    return kept, stats


def main():
    args = parse_args()

    # Run QC
    kept, stats = qc_filter(args.input, args.min_len, args.max_ambig)

    # Write clean FASTA
    SeqIO.write(kept, args.output, "fasta")

    # Write stats CSV
    with open(args.stats, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "input", "kept",
            "removed_short", "removed_ambig", "removed_frame",
            "removed_stop", "removed_dup"
        ])
        writer.writerow([
            args.gene, stats["input"], stats["kept"],
            stats["removed_short"], stats["removed_ambig"],
            stats["removed_frame"], stats["removed_stop"],
            stats["removed_dup"]
        ])

    # Print summary
    print(f"[QC] {args.gene}: {stats['input']} input → {stats['kept']} kept")
    print(f"  Removed: short={stats['removed_short']}, ambig={stats['removed_ambig']}, "
          f"frame={stats['removed_frame']}, stop={stats['removed_stop']}, "
          f"dup={stats['removed_dup']}")


if __name__ == "__main__":
    main()
