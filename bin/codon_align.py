#!/usr/bin/env python3
"""
Protein-guided codon alignment.

Steps:
  1. Translate nucleotide → amino acid
  2. Align protein sequences with MAFFT (--auto)
  3. Back-translate to codon-aware nucleotide alignment

This ensures the alignment respects codon boundaries — critical for
downstream dN/dS selection analysis.

Usage:
  codon_align.py --input gene.clean.fasta --output gene.codon.fasta
"""

import argparse
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description="Protein-guided codon alignment")
    parser.add_argument("--input", required=True, help="Input nucleotide FASTA")
    parser.add_argument("--output", required=True, help="Output codon alignment FASTA")
    parser.add_argument("--gene", default="gene", help="Gene name for logging")
    parser.add_argument("--threads", type=int, default=2, help="MAFFT threads")
    return parser.parse_args()


def translate_sequences(nt_records):
    """Translate nucleotide sequences to amino acids (table 11 = bacterial)."""
    aa_records = []
    skipped = 0
    for rec in nt_records:
        seq = str(rec.seq).upper()
        # Ensure length is divisible by 3
        if len(seq) % 3 != 0:
            skipped += 1
            continue
        try:
            aa_seq = Seq(seq).translate(table=11, to_stop=False)
            # Remove trailing stop codon if present
            aa_str = str(aa_seq).rstrip("*")
            aa_records.append(
                SeqRecord(Seq(aa_str), id=rec.id, description="")
            )
        except Exception as e:
            print(f"[warn] Could not translate {rec.id}: {e}", file=sys.stderr)
            skipped += 1
    if skipped:
        print(f"[align] Skipped {skipped} sequences during translation")
    return aa_records


def run_mafft(input_fasta, output_fasta, threads=2):
    """Run MAFFT alignment on protein sequences."""
    cmd = [
        "mafft",
        "--auto",
        "--thread", str(threads),
        "--quiet",
        input_fasta
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[error] MAFFT failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    with open(output_fasta, "w") as f:
        f.write(result.stdout)

    return output_fasta


def backtranslate(aa_alignment, nt_records):
    """Back-translate protein alignment to codon alignment."""
    # Build lookup: id → nucleotide sequence
    id2nt = {rec.id: str(rec.seq).upper() for rec in nt_records}

    codon_records = []
    for aa_rec in aa_alignment:
        if aa_rec.id not in id2nt:
            print(f"[warn] No NT sequence for {aa_rec.id}", file=sys.stderr)
            continue

        nt_seq = id2nt[aa_rec.id]
        codons = [nt_seq[i:i+3] for i in range(0, len(nt_seq), 3)]

        # Map aligned AA to codons
        codon_aln = []
        codon_idx = 0
        for aa in str(aa_rec.seq):
            if aa == "-":
                codon_aln.append("---")
            else:
                if codon_idx < len(codons):
                    codon_aln.append(codons[codon_idx])
                else:
                    codon_aln.append("---")
                codon_idx += 1

        codon_records.append(
            SeqRecord(Seq("".join(codon_aln)), id=aa_rec.id, description="")
        )

    return codon_records


def main():
    args = parse_args()

    # Step 1: Read input nucleotide sequences
    nt_records = list(SeqIO.parse(args.input, "fasta"))
    print(f"[{args.gene}] Input: {len(nt_records)} sequences")

    if len(nt_records) < 2:
        print(f"[{args.gene}] Too few sequences for alignment, copying as-is")
        SeqIO.write(nt_records, args.output, "fasta")
        return

    # Step 2: Translate to amino acids (using bacterial genetic code, table 11)
    aa_records = translate_sequences(nt_records)
    print(f"[{args.gene}] Translated: {len(aa_records)} AA sequences")

    # Step 3: Write AA to temp file and align with MAFFT
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
        aa_tmp = f.name
        SeqIO.write(aa_records, f, "fasta")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".aln.faa", delete=False) as f:
        aa_aln_tmp = f.name

    run_mafft(aa_tmp, aa_aln_tmp, threads=args.threads)
    aa_alignment = list(SeqIO.parse(aa_aln_tmp, "fasta"))
    print(f"[{args.gene}] MAFFT aligned: {len(aa_alignment)} sequences, "
          f"length {len(str(aa_alignment[0].seq)) if aa_alignment else 0} aa")

    # Step 4: Back-translate to codon alignment
    codon_records = backtranslate(aa_alignment, nt_records)

    # Step 5: Write codon alignment
    SeqIO.write(codon_records, args.output, "fasta")
    aln_len = len(str(codon_records[0].seq)) if codon_records else 0
    print(f"[{args.gene}] Codon alignment: {len(codon_records)} sequences, {aln_len} nt")


if __name__ == "__main__":
    main()
