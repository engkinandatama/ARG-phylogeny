#!/usr/bin/env python3
"""
Concatenate multiple codon alignments by common taxa (intersection).

Usage:
  concat_alignments.py --input-dir . --output concat.fasta --pattern "*.codon.fasta"
"""

import argparse
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description="Concatenate alignments by common taxa")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pattern", default="*.codon.fasta")
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    aln_files = sorted(input_dir.glob(args.pattern))

    if not aln_files:
        print("[concat] No alignment files found!")
        open(args.output, "w").close()
        return

    # Find common taxa across all alignments
    taxa_sets = []
    for aln_file in aln_files:
        taxa = set(r.id for r in SeqIO.parse(str(aln_file), "fasta"))
        taxa_sets.append(taxa)

    common_taxa = set.intersection(*taxa_sets) if taxa_sets else set()
    print(f"[concat] {len(aln_files)} alignments, {len(common_taxa)} common taxa")

    if not common_taxa:
        print("[concat] WARNING: No common taxa found across alignments!")
        open(args.output, "w").close()
        return

    # Concatenate sequences for common taxa
    concat_seqs = defaultdict(list)
    for aln_file in aln_files:
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(str(aln_file), "fasta")}
        for taxon in sorted(common_taxa):
            concat_seqs[taxon].append(seqs.get(taxon, ""))

    # Write concatenated alignment
    records = []
    for taxon in sorted(common_taxa):
        concat_seq = "".join(concat_seqs[taxon])
        records.append(SeqRecord(Seq(concat_seq), id=taxon, description=""))

    SeqIO.write(records, args.output, "fasta")
    total_len = len(str(records[0].seq)) if records else 0
    print(f"[concat] Output: {len(records)} taxa, {total_len} nt total length")


if __name__ == "__main__":
    main()
