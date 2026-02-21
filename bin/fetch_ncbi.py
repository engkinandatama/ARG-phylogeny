#!/usr/bin/env python3
"""
Fetch CDS sequences from NCBI GenBank using Entrez queries.

Outputs:
  - FASTA file with rich headers (accession|gene|product|protein_id|country|date)
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
import sys
import time

from Bio import Entrez, SeqIO


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
    return parser.parse_args()


def fetch_ids(query, max_records, retries=3):
    """Search NCBI and return list of IDs."""
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_records)
            result = Entrez.read(handle)
            handle.close()
            ids = result["IdList"]
            print(f"[fetch] Found {len(ids)} records (query: {query[:80]}...)")
            return ids
        except Exception as e:
            print(f"[fetch] esearch attempt {attempt+1} failed: {e}", file=sys.stderr)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Failed to search NCBI after {retries} attempts")


def fetch_records(ids, batch_size=100, retries=3):
    """Fetch GenBank records in batches."""
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
                print(f"[fetch] Batch {start//batch_size + 1}: {len(records)} records")
                time.sleep(0.5)  # Be nice to NCBI
                break
            except Exception as e:
                print(f"[fetch] efetch attempt {attempt+1} failed: {e}", file=sys.stderr)
                time.sleep(5 * (attempt + 1))
        else:
            print(f"[fetch] WARNING: Skipping batch starting at {start}", file=sys.stderr)
    return all_records


def extract_cds(records, gene_name, min_cds_len=90):
    """Extract CDS features from GenBank records."""
    rows = []
    for rec in records:
        cds_features = [f for f in rec.features if f.type == "CDS"]
        for cds in cds_features:
            seq = str(cds.extract(rec.seq)).upper()

            # Skip short sequences
            if len(seq) < min_cds_len:
                continue

            # Extract metadata from qualifiers
            gene = cds.qualifiers.get("gene", [""])[0]
            product = cds.qualifiers.get("product", [""])[0]
            protein_id = cds.qualifiers.get("protein_id", [""])[0]

            # Extract metadata from record annotations
            source_features = [f for f in rec.features if f.type == "source"]
            country = ""
            isolation_source = ""
            if source_features:
                country = source_features[0].qualifiers.get("country", [""])[0]
                isolation_source = source_features[0].qualifiers.get("isolation_source", [""])[0]

            # Collection date from annotations or source feature
            date = rec.annotations.get("date", "")
            if source_features and not date:
                date = source_features[0].qualifiers.get("collection_date", [""])[0]

            # Build header with pipe-separated metadata
            header = (
                f"{rec.id}|gene={gene}|product={product}"
                f"|protein_id={protein_id}|country={country}|date={date}"
            )

            rows.append({
                "header": header,
                "seq": seq,
                "acc": rec.id,
                "gene": gene,
                "product": product,
                "protein_id": protein_id,
                "country": country,
                "isolation_source": isolation_source,
                "date": date,
                "length_nt": len(seq),
            })

    print(f"[{gene_name}] Extracted {len(rows)} CDS sequences")
    return rows


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

    # Step 1: Search NCBI
    ids = fetch_ids(args.query, args.max_records)
    if not ids:
        print(f"[{args.gene}] No records found. Creating empty outputs.")
        open(args.out_fasta, "w").close()
        with open(args.out_meta, "w") as f:
            f.write("acc,gene,product,protein_id,country,isolation_source,date,length_nt\n")
        return

    # Step 2: Fetch GenBank records
    records = fetch_records(ids, batch_size=args.batch_size)

    # Step 3: Extract CDS
    rows = extract_cds(records, args.gene, min_cds_len=args.min_cds_len)

    # Step 4: Write outputs
    write_outputs(rows, args.out_fasta, args.out_meta)


if __name__ == "__main__":
    main()
