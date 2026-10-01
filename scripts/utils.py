"""
utils.py — shared helpers used across the pipeline.

Import from other scripts with:
    from utils import parse_header, load_metadata, decimal_year, dec_to_dt
"""

import csv
import math
from datetime import datetime, timedelta


# ── FASTA header format: >accession|YYYY-MM-DD|Country_name ─────────────────

def parse_header(raw_header: str) -> dict:
    """
    Parse the pipeline's FASTA header format.

    '>MH708882|2018-06-15|USA' → {'accession': 'MH708882',
                                   'date_str':  '2018-06-15',
                                   'country':   'USA'}
    Missing fields are returned as 'unknown'.
    The leading '>' is stripped automatically.
    """
    parts = raw_header.lstrip(">").split("|")
    return {
        "accession": parts[0] if len(parts) > 0 else "unknown",
        "date_str":  parts[1] if len(parts) > 1 else "unknown",
        "country":   parts[2] if len(parts) > 2 else "unknown",
    }


def parse_date_str(date_str: str) -> datetime | None:
    """
    Parse an ISO-8601 date string to a datetime.
    Returns None on failure.
    Handles 'YYYY-MM-DD', 'YYYY-MM', 'YYYY'.
    """
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None


# ── metadata CSV (written by fetch_sequences.py) ─────────────────────────────

def load_metadata(csv_path: str) -> dict:
    """
    Read data/metadata.csv into a dict keyed by accession.

    Returns:
        { accession: {'date': datetime, 'year': int, 'country': str,
                      'vp1_length': int} }
    """
    meta = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            date_obj = parse_date_str(row.get("date", ""))
            if date_obj is None:
                continue
            meta[row["accession"]] = {
                "date":       date_obj,
                "year":       int(row.get("year", date_obj.year)),
                "country":    row.get("country", "unknown"),
                "vp1_length": int(row.get("vp1_length", 0)),
            }
    return meta


# ── date arithmetic ──────────────────────────────────────────────────────────

def decimal_year(dt: datetime) -> float:
    """
    Convert a datetime to a decimal year.

    datetime(2014, 7, 2) → 2014.501 (roughly)
    """
    year_start = datetime(dt.year, 1, 1)
    year_end   = datetime(dt.year + 1, 1, 1)
    fraction   = (dt - year_start) / (year_end - year_start)
    return dt.year + fraction


def dec_to_dt(dy: float) -> datetime:
    """
    Inverse of decimal_year — convert a decimal year back to a datetime.

    2014.5 → datetime(2014, 7, 2)
    """
    yr    = int(dy)
    frac  = dy - yr
    days  = (datetime(yr + 1, 1, 1) - datetime(yr, 1, 1)).days
    return datetime(yr, 1, 1) + timedelta(days=frac * days)


# ── phylogenetic tree helpers ─────────────────────────────────────────────────

def accession_from_tip(tip_name: str) -> str:
    """
    Strip date and country from a tip name to recover the accession.

    'MH708882|2018-06-15|USA' → 'MH708882'
    """
    return tip_name.split("|")[0]


def root_to_tip_distances(tree) -> dict:
    """
    Return { tip_name: root_to_tip_distance } for a rooted BioPython Tree.
    """
    root    = tree.root
    tip_map = {}
    for tip in tree.get_terminals():
        tip_map[tip.name] = tree.distance(root, tip)
    return tip_map
