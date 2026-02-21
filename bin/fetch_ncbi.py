#!/usr/bin/env python3
"""
Fetch CDS sequences from NCBI GenBank using Entrez queries.

Two-pass strategy for handling both small records and WGS contigs:
  1. GB records → metadata (country, date) + sequences for small records
  2. fasta_cds_na fallback → CDS sequences for WGS/large records

Filters:
  - Excludes computationally predicted (PREDICTED, MODEL REFSEQ)
  - Excludes unverified and synthetic constructs
  - Excludes metagenomic/environmental samples
  - Minimum CDS length filter

Outputs:
  - FASTA file with rich headers
  - CSV metadata file

Usage:
  fetch_ncbi.py --gene blaCTX-M \
                --query "(Klebsiella pneumoniae[Organism]) AND blaCTX-M[Gene]..." \
                --max-records 200 \
                --email you@email.com \
                --out-fasta blaCTX-M.fasta \
                --out-meta blaCTX-M_meta.csv
"""

import argparse
import csv
import re
import sys
import time

from Bio import Entrez, SeqIO


# ── Exclusion filters for computationally predicted sequences ─────────
EXCLUDE_TERMS = [
    'NOT "PREDICTED"[Title]',
    'NOT "UNVERIFIED"[Title]',
    'NOT "synthetic construct"[Organism]',
    'NOT "synthetic"[Title]',
    'NOT metagenome[Organism]',
    'NOT "environmental sample"[Isolation Source]',
]


def parse_args():
    parser = argparse.ArgumentParser(description="Fetch CDS from NCBI GenBank")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--query", required=True, help="NCBI Entrez query")
    parser.add_argument("--max-records", type=int, default=200, help="Max records to fetch")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--out-fasta", required=True, help="Output FASTA path")
    parser.add_argument("--out-meta", required=True, help="Output metadata CSV path")
    parser.add_argument("--min-cds-len", type=int, default=90, help="Min CDS length in nt")
    parser.add_argument("--batch-size", type=int, default=100, help="Records per Entrez batch")
    parser.add_argument("--no-filter", action="store_true",
                        help="Disable exclusion of predicted/synthetic sequences")
    return parser.parse_args()


def build_filtered_query(base_query, apply_filter=True):
    """Add exclusion filters to query to remove predicted/synthetic records."""
    if not apply_filter:
        return base_query
    parts = [base_query] + EXCLUDE_TERMS
    return " ".join(parts)


def fetch_ids(query, max_records, retries=3):
    """Search NCBI and return list of IDs."""
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_records)
            result = Entrez.read(handle)
            handle.close()
            ids = result["IdList"]
            total = int(result.get("Count", len(ids)))
            print(f"[fetch] Found {total} total, fetching {len(ids)} records")
            return ids
        except Exception as e:
            print(f"[fetch] esearch attempt {attempt+1} failed: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Failed to search NCBI after {retries} attempts")


def fetch_records_gb(ids, batch_size=100, retries=3):
    """Fetch GenBank records in batches (Pass 1: metadata + inline sequences)."""
    all_records = []
    for start in range(0, len(ids), batch_size):
        batch = ids[start:start + batch_size]
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(batch),
                    rettype="gb",
                    retmode="text"
                )
                records = list(SeqIO.parse(handle, "gb"))
                handle.close()
                all_records.extend(records)
                print(f"[fetch] GB batch {start//batch_size + 1}: {len(records)} records")
                time.sleep(0.5)  # Be nice to NCBI
                break
            except Exception as e:
                print(f"[fetch] efetch attempt {attempt+1} failed: {e}", file=sys.stderr)
                time.sleep(5 * (attempt + 1))
        else:
            print(f"[fetch] WARNING: Skipping batch starting at {start}", file=sys.stderr)
    return all_records


def gene_matches(gene_qualifier, target_gene):
    """Check if gene qualifier matches target gene (case-insensitive, variant-aware).

    Handles cases like:
      - blaCTX-M-15 matches blaCTX-M
      - blaKPC-2 matches blaKPC
      - oqxA1 matches oqxA
      - rpoB matches rpoB
    """
    if not gene_qualifier:
        return False
    g = gene_qualifier.lower().strip()
    t = target_gene.lower().strip()
    # Exact match or target is a prefix (e.g., blaCTX-M matches blaCTX-M-15)
    return g == t or g.startswith(t)


def extract_metadata(rec):
    """Extract metadata from a GenBank record (country, date, isolation source)."""
    source_features = [f for f in rec.features if f.type == "source"]
    country = ""
    isolation_source = ""
    if source_features:
        country = source_features[0].qualifiers.get("country", [""])[0]
        isolation_source = source_features[0].qualifiers.get("isolation_source", [""])[0]

    date = rec.annotations.get("date", "")
    if source_features and not date:
        date = source_features[0].qualifiers.get("collection_date", [""])[0]

    return {
        "country": country,
        "isolation_source": isolation_source,
        "date": date,
    }


def extract_from_gb(records, gene_name, min_cds_len=90):
    """Pass 1: Extract CDS from GB records with inline sequences.

    Returns:
        rows: list of dicts with extracted data
        undefined_ids: list of record IDs that had undefined sequences
        metadata_by_acc: dict of metadata keyed by accession (for Pass 2)
    """
    rows = []
    undefined_ids = []
    metadata_by_acc = {}

    for rec in records:
        # Always extract metadata (works even with undefined sequences)
        meta = extract_metadata(rec)
        metadata_by_acc[rec.id] = meta

        has_undefined = False
        for cds in [f for f in rec.features if f.type == "CDS"]:
            # Filter by gene name first
            gene = cds.qualifiers.get("gene", [""])[0]
            if not gene_matches(gene, gene_name):
                continue

            # Try to extract sequence
            try:
                seq = str(cds.extract(rec.seq)).upper()
            except Exception:
                has_undefined = True
                continue

            if len(seq) < min_cds_len:
                continue

            product = cds.qualifiers.get("product", [""])[0]
            protein_id = cds.qualifiers.get("protein_id", [""])[0]

            header = (
                f"{rec.id}|gene={gene}|product={product}"
                f"|protein_id={protein_id}|country={meta['country']}|date={meta['date']}"
            )

            rows.append({
                "header": header,
                "seq": seq,
                "acc": rec.id,
                "gene": gene,
                "product": product,
                "protein_id": protein_id,
                "country": meta["country"],
                "isolation_source": meta["isolation_source"],
                "date": meta["date"],
                "length_nt": len(seq),
            })

        if has_undefined:
            undefined_ids.append(rec.id)

    return rows, undefined_ids, metadata_by_acc


def parse_fasta_field(description, field):
    """Parse [field=value] from NCBI fasta_cds_na FASTA description."""
    match = re.search(rf'\[{field}=([^\]]+)\]', description)
    return match.group(1) if match else ""


def parse_parent_accession(fasta_id):
    """Extract parent accession from fasta_cds_na record ID.

    Example: lcl|NZ_CP012345.1_cds_WP_004152123.1_1234 → NZ_CP012345.1
    """
    match = re.match(r'lcl\|(\S+?)_cds_', fasta_id)
    return match.group(1) if match else fasta_id


def fetch_cds_na_fallback(ids, gene_name, metadata_by_acc, min_cds_len=90,
                          batch_size=20, retries=3):
    """Pass 2: Fetch CDS nucleotide FASTA for records with undefined sequences.

    Uses rettype=fasta_cds_na which tells NCBI to extract CDS server-side,
    avoiding the need to download entire WGS contigs.
    """
    rows = []
    total_cds_checked = 0

    # Smaller batches for fasta_cds_na (returns many entries per record)
    for start in range(0, len(ids), batch_size):
        batch = ids[start:start + batch_size]
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(batch),
                    rettype="fasta_cds_na",
                    retmode="text"
                )

                for rec in SeqIO.parse(handle, "fasta"):
                    total_cds_checked += 1

                    # Filter by gene name from FASTA header
                    gene = parse_fasta_field(rec.description, "gene")
                    if not gene_matches(gene, gene_name):
                        continue

                    seq = str(rec.seq).upper()
                    if len(seq) < min_cds_len:
                        continue

                    # Parse metadata from FASTA header
                    protein_id = parse_fasta_field(rec.description, "protein_id")
                    product = parse_fasta_field(rec.description, "protein")
                    acc = parse_parent_accession(rec.id)

                    # Enrich with metadata from GB records (Pass 1)
                    meta = metadata_by_acc.get(acc, {})
                    country = meta.get("country", "")
                    isolation_source = meta.get("isolation_source", "")
                    date = meta.get("date", "")

                    header = (
                        f"{acc}|gene={gene}|product={product}"
                        f"|protein_id={protein_id}|country={country}|date={date}"
                    )

                    rows.append({
                        "header": header,
                        "seq": seq,
                        "acc": acc,
                        "gene": gene,
                        "product": product,
                        "protein_id": protein_id,
                        "country": country,
                        "isolation_source": isolation_source,
                        "date": date,
                        "length_nt": len(seq),
                    })

                handle.close()
                print(f"[fetch] CDS-NA batch {start//batch_size + 1}: "
                      f"checked {total_cds_checked} CDSs, found {len(rows)} matches")
                time.sleep(0.5)
                break
            except Exception as e:
                print(f"[fetch] fasta_cds_na attempt {attempt+1} failed: {e}",
                      file=sys.stderr)
                time.sleep(5 * (attempt + 1))
        else:
            print(f"[fetch] WARNING: Skipping CDS-NA batch at {start}", file=sys.stderr)

    return rows


def deduplicate_rows(rows):
    """Remove duplicate sequences (same accession + same sequence)."""
    seen = set()
    unique = []
    for row in rows:
        key = (row["acc"], row["seq"])
        if key not in seen:
            seen.add(key)
            unique.append(row)
    return unique


def write_outputs(rows, out_fasta, out_meta):
    """Write FASTA and metadata CSV files."""
    # Write FASTA
    with open(out_fasta, "w") as f:
        for row in rows:
            f.write(f">{row['header']}\n{row['seq']}\n")

    # Write metadata CSV
    fieldnames = [
        "acc", "gene", "product", "protein_id",
        "country", "isolation_source", "date", "length_nt"
    ]
    with open(out_meta, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"[output] FASTA: {out_fasta} ({len(rows)} sequences)")
    print(f"[output] Meta:  {out_meta}")


def main():
    args = parse_args()
    Entrez.email = args.email

    # Build filtered query (exclude predicted/synthetic)
    query = build_filtered_query(args.query, apply_filter=not args.no_filter)
    print(f"[{args.gene}] Query: {query[:120]}...")

    # Step 1: Search NCBI
    ids = fetch_ids(query, args.max_records)
    if not ids:
        print(f"[{args.gene}] No records found. Creating empty outputs.")
        open(args.out_fasta, "w").close()
        with open(args.out_meta, "w") as f:
            f.write("acc,gene,product,protein_id,country,isolation_source,date,length_nt\n")
        return

    # Step 2: Pass 1 — Fetch GB records (metadata + inline sequences)
    print(f"[{args.gene}] Pass 1: Fetching GenBank records...")
    gb_records = fetch_records_gb(ids, batch_size=args.batch_size)

    # Step 3: Extract CDS from GB records
    rows, undefined_ids, metadata_by_acc = extract_from_gb(
        gb_records, args.gene, min_cds_len=args.min_cds_len
    )
    print(f"[{args.gene}] Pass 1 result: {len(rows)} sequences extracted, "
          f"{len(undefined_ids)} records had undefined sequences")

    # Step 4: Pass 2 — Fallback for records with undefined sequences (WGS contigs)
    if undefined_ids:
        print(f"[{args.gene}] Pass 2: Re-fetching {len(undefined_ids)} WGS records "
              f"via CDS extraction...")
        fallback_rows = fetch_cds_na_fallback(
            undefined_ids, args.gene, metadata_by_acc,
            min_cds_len=args.min_cds_len
        )
        print(f"[{args.gene}] Pass 2 result: {len(fallback_rows)} additional sequences")
        rows.extend(fallback_rows)

    # Step 5: Deduplicate
    rows = deduplicate_rows(rows)
    print(f"[{args.gene}] Final: {len(rows)} unique sequences")

    # Step 6: Write outputs
    write_outputs(rows, args.out_fasta, args.out_meta)


if __name__ == "__main__":
    main()
