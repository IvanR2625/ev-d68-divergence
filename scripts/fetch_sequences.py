"""
Step 1 — Fetch EV-D68 VP1 sequences from NCBI.

Queries NCBI Nucleotide for complete Enterovirus D68 genomes,
extracts the VP1 coding sequence from each record, parses the
collection date, and writes:
  data/ev_d68_vp1_raw.fasta   — one VP1 sequence per sample
  data/metadata.csv            — accession, date, country, length

Run this first. Expects an internet connection.
"""

import os
import re
import sys
import time
import csv
import random
from datetime import datetime
from collections import defaultdict

from Bio import Entrez, SeqIO

# ── configuration ────────────────────────────────────────────────────────────
# Replace with your real email — NCBI requires it for Entrez API access.
ENTREZ_EMAIL = "your.email@example.com"

# Years to cover. Adjust the end year as needed.
YEAR_START = 2014
YEAR_END   = 2024

# Maximum sequences to keep per year (random sample if more are found).
# 15 × 11 years = 165 sequences — a manageable MAFFT input.
MAX_PER_YEAR = 15

# Output paths (relative to the project root)
OUT_FASTA    = os.path.join("data", "ev_d68_vp1_raw.fasta")
OUT_METADATA = os.path.join("data", "metadata.csv")
# ─────────────────────────────────────────────────────────────────────────────


def parse_collection_date(gb_record):
    """
    Extract and normalise the collection_date qualifier from the SOURCE
    feature of a GenBank record.

    Returns a datetime, or None if the date cannot be parsed.
    GenBank uses many formats: '2014-09-15', '2014-09', '2014',
    'Sep-2014', '15-Sep-2014'.
    """
    for feature in gb_record.features:
        if feature.type == "source":
            if "collection_date" in feature.qualifiers:
                raw = feature.qualifiers["collection_date"][0]
                # Try formats from most to least specific.
                for fmt in ("%Y-%m-%d", "%d-%b-%Y", "%b-%Y", "%Y-%m", "%Y"):
                    try:
                        return datetime.strptime(raw, fmt)
                    except ValueError:
                        continue
    return None


def parse_country(gb_record):
    """Return the country string from the source feature, or 'unknown'."""
    for feature in gb_record.features:
        if feature.type == "source":
            if "country" in feature.qualifiers:
                # Strip sub-region (e.g. 'USA: California' → 'USA')
                return feature.qualifiers["country"][0].split(":")[0].strip()
    return "unknown"


def extract_vp1(gb_record):
    """
    Find the VP1 CDS feature and return its nucleotide sequence as a
    string, or None if the annotation is absent.

    EV-D68 genomes on GenBank annotate VP1 with gene='VP1' and
    feature type CDS.
    """
    for feature in gb_record.features:
        if feature.type != "CDS":
            continue
        gene_name = feature.qualifiers.get("gene", [""])[0].upper()
        product   = feature.qualifiers.get("product", [""])[0].upper()
        if "VP1" in gene_name or "VP1" in product:
            # Extract the sub-sequence from the full record.
            vp1_seq = feature.extract(gb_record.seq)
            return str(vp1_seq)
    return None


def fetch_accessions_for_year(year):
    """
    Return a list of NCBI accession numbers for complete EV-D68 genomes
    collected in `year`, up to a generous ceiling of 500.
    """
    query = (
        f'"Enterovirus D68"[Organism] '
        f'AND "{year}/01/01"[PDAT] : "{year}/12/31"[PDAT] '
        f'AND "complete genome"[Title]'
    )
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=500)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


def fetch_records(id_list, batch_size=50):
    """
    Fetch GenBank records for a list of GI/accession IDs in batches.
    Yields one SeqRecord at a time.
    """
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i : i + batch_size]
        ids   = ",".join(batch)
        handle = Entrez.efetch(
            db="nucleotide", id=ids, rettype="gb", retmode="text"
        )
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()
        for rec in records:
            yield rec
        # Be polite to NCBI — max ~3 requests/second without an API key.
        time.sleep(0.4)


def safe_header(accession, date_obj, country):
    """Build a FASTA header that encodes accession, ISO date, and country."""
    date_str = date_obj.strftime("%Y-%m-%d") if date_obj else "unknown"
    # Strip spaces from country so the header parses cleanly later.
    country_clean = country.replace(" ", "_")
    return f"{accession}|{date_str}|{country_clean}"


def main():
    if ENTREZ_EMAIL == "your.email@example.com":
        sys.exit(
            "ERROR: Set ENTREZ_EMAIL at the top of this script to your "
            "real email address before running. NCBI requires it."
        )

    Entrez.email = ENTREZ_EMAIL

    os.makedirs("data", exist_ok=True)

    records_by_year = defaultdict(list)
    total_found = 0

    print(f"Searching NCBI for EV-D68 complete genomes, {YEAR_START}–{YEAR_END} …")

    for year in range(YEAR_START, YEAR_END + 1):
        print(f"  {year}: ", end="", flush=True)
        id_list = fetch_accessions_for_year(year)
        print(f"{len(id_list)} accessions found", end="")

        if not id_list:
            print()
            continue

        # Sample down to MAX_PER_YEAR if needed (reproducible with seed).
        random.seed(year)
        if len(id_list) > MAX_PER_YEAR:
            id_list = random.sample(id_list, MAX_PER_YEAR)
            print(f", sampling {MAX_PER_YEAR}")
        else:
            print(f", keeping all {len(id_list)}")

        for gb_rec in fetch_records(id_list):
            date_obj = parse_collection_date(gb_rec)
            country  = parse_country(gb_rec)
            vp1_seq  = extract_vp1(gb_rec)

            if vp1_seq is None:
                continue          # skip records with no VP1 annotation
            if date_obj is None:
                continue          # skip records with no collection date

            records_by_year[year].append({
                "accession": gb_rec.id,
                "date":      date_obj,
                "country":   country,
                "vp1_seq":   vp1_seq,
            })
            total_found += 1

        time.sleep(0.5)   # extra pause between years

    print(f"\n{total_found} sequences with VP1 annotation and collection date.")

    if total_found == 0:
        sys.exit("No usable sequences found. Check your search terms or network.")

    # ── write FASTA ──────────────────────────────────────────────────────────
    fasta_rows = []
    meta_rows  = []

    for year, entries in sorted(records_by_year.items()):
        for e in entries:
            header = safe_header(e["accession"], e["date"], e["country"])
            fasta_rows.append(f">{header}")
            fasta_rows.append(e["vp1_seq"])
            meta_rows.append({
                "accession": e["accession"],
                "date":      e["date"].strftime("%Y-%m-%d"),
                "year":      e["date"].year,
                "country":   e["country"],
                "vp1_length": len(e["vp1_seq"]),
            })

    with open(OUT_FASTA, "w") as fh:
        fh.write("\n".join(fasta_rows) + "\n")

    with open(OUT_METADATA, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["accession", "date", "year", "country", "vp1_length"]
        )
        writer.writeheader()
        writer.writerows(meta_rows)

    print(f"Wrote {OUT_FASTA}  ({total_found} sequences)")
    print(f"Wrote {OUT_METADATA}")


if __name__ == "__main__":
    main()
