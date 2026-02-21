#!/usr/bin/env python3
"""
Generate synthetic test data for pipeline validation.
Creates small FASTA files + metadata CSVs that mimic real NCBI data.

Usage:
  python generate_test_data.py --outdir test/test_data
"""

import argparse
import csv
import random
import string
from pathlib import Path

# Realistic coding sequences (divisible by 3, no internal stops in table 11)
# These are synthetic but structurally valid
CODON_TABLE = [
    "ATG", "GCT", "GCC", "GCA", "AAA", "AAG", "GAT", "GAC",
    "TGT", "TGC", "GAA", "GAG", "TTT", "TTC", "GGT", "GGC",
    "CAT", "CAC", "ATT", "ATC", "CTG", "CTA", "CTT", "CTC",
    "ATG", "AAT", "AAC", "CCT", "CCC", "CAA", "CAG", "CGT",
    "CGC", "TCT", "TCC", "ACT", "ACC", "GTT", "GTC", "TGG",
    "TAT", "TAC", "GCG", "GGA", "GGG", "CGA", "CGG", "AGT",
    "AGC", "AGA", "AGG", "TTA", "TTG", "GTA", "GTG", "CCA",
    "CCG", "ACA", "ACG", "TCA", "TCG",
]

COUNTRIES = [
    "China", "India", "USA", "Brazil", "Thailand",
    "South Africa", "Italy", "Japan", "Nigeria", "Egypt"
]

YEARS = list(range(2015, 2025))


def random_cds(length_codons=280):
    """Generate a random CDS (starts with ATG, ends with stop, no internal stops)."""
    seq = ["ATG"]  # start codon
    for _ in range(length_codons - 2):
        seq.append(random.choice(CODON_TABLE))
    seq.append("TAA")  # stop codon
    return "".join(seq)


def introduce_variants(base_seq, n_variants, mutation_rate=0.02):
    """Create n variants from a base sequence with small mutations."""
    variants = []
    seq_list = list(base_seq)
    for i in range(n_variants):
        var = seq_list.copy()
        n_mutations = max(1, int(len(var) * mutation_rate * (i + 1) / n_variants))
        for _ in range(n_mutations):
            pos = random.randint(0, len(var) - 1)
            var[pos] = random.choice("ACGT")
        # Ensure still valid: starts ATG, divisible by 3
        var[0], var[1], var[2] = "A", "T", "G"
        seq = "".join(var)
        # Trim to codon boundary
        seq = seq[:len(seq) - len(seq) % 3]
        variants.append(seq)
    return variants


def generate_gene_data(gene_name, n_seqs, outdir):
    """Generate FASTA + metadata CSV for a gene."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Base sequence length depends on gene type
    lengths = {
        "blaCTX-M": 290, "blaSHV": 286, "blaTEM": 286,
        "blaKPC": 293, "blaNDM": 270, "blaOXA-48": 265,
        "oqxA": 400, "qnrB": 214, "sul1": 280,
        "tetA": 400, "mcr-1": 541,
        "rpoB": 1100, "gyrA": 880, "gapA": 330,
        "rpoD": 560, "infB": 270,
    }
    n_codons = lengths.get(gene_name, 300)

    # Generate base sequence + variants
    base = random_cds(n_codons)
    variants = [base] + introduce_variants(base, n_seqs - 1)

    # Write FASTA
    fasta_path = outdir / f"{gene_name}.fasta"
    with open(fasta_path, "w") as f:
        for i, seq in enumerate(variants):
            acc = f"TEST{i+1:03d}"
            f.write(f">{acc}|gene={gene_name}|country={random.choice(COUNTRIES)}\n")
            # Wrap at 70 chars
            for j in range(0, len(seq), 70):
                f.write(seq[j:j+70] + "\n")

    # Write metadata CSV
    meta_path = outdir / f"{gene_name}_meta.csv"
    with open(meta_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["acc", "gene", "product", "protein_id", "country", "date", "length_nt"])
        for i, seq in enumerate(variants):
            acc = f"TEST{i+1:03d}"
            country = random.choice(COUNTRIES)
            year = random.choice(YEARS)
            w.writerow([
                acc, gene_name, f"{gene_name} protein", f"PID_{acc}",
                country, f"{year}-{random.randint(1,12):02d}-{random.randint(1,28):02d}",
                len(seq)
            ])

    print(f"[test] {gene_name}: {n_seqs} sequences ({len(base)} nt) → {fasta_path}")
    return fasta_path, meta_path


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic test data")
    parser.add_argument("--outdir", default="test/test_data")
    parser.add_argument("--n-seqs", type=int, default=8)
    args = parser.parse_args()

    # Mini test genes (matching test_mini.yaml)
    test_genes = ["blaCTX-M", "blaKPC", "rpoB", "gyrA"]

    print(f"[test] Generating synthetic data for {len(test_genes)} genes...")
    for gene in test_genes:
        generate_gene_data(gene, args.n_seqs, args.outdir)

    print(f"\n[test] Done! Files saved to: {args.outdir}/")
    print(f"[test] Run pipeline with: nextflow run main.nf -params-file conf/test_mini.yaml")


if __name__ == "__main__":
    main()
