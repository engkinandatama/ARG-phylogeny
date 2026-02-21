#!/usr/bin/env python3
"""
CARD/ResFinder variant naming via BLAST.

Placeholder: When CARD database is available, runs BLAST against it.
Currently generates a name based on gene family + accession.

Usage:
  card_lookup.py --input gene.fasta --gene blaCTX-M --output variants.csv
"""

import argparse
import csv

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="CARD variant naming")
    parser.add_argument("--input", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--card-db", default=None, help="Path to CARD BLAST DB")
    return parser.parse_args()


def main():
    args = parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))

    results = []
    for rec in records:
        # TODO: When CARD DB is available, run BLAST:
        # blastn -query input.fasta -db card_db -outfmt 6 -max_target_seqs 1
        # Parse best hit for variant name

        results.append({
            "accession": rec.id,
            "gene": args.gene,
            "variant": f"{args.gene}-like",  # Placeholder
            "identity": 100.0,
        })

    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "gene", "variant", "identity"])
        w.writeheader()
        w.writerows(results)

    print(f"[card] {args.gene}: {len(results)} sequences processed")


if __name__ == "__main__":
    main()
