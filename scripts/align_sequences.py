"""
Step 2 — Align VP1 sequences with MAFFT.

Reads  data/ev_d68_vp1_raw.fasta
Writes data/ev_d68_vp1_aligned.fasta

Requires MAFFT to be installed and on PATH.
Download: https://mafft.cbrc.jp/alignment/software/

On Windows:
  1. Download the Windows installer from the link above.
  2. Install it (default path: C:\Program Files\MAFFT\).
  3. Add the install folder to your system PATH, or set MAFFT_EXE below.

MAFFT flags used:
  --auto        let MAFFT choose the best algorithm for the input size
  --thread -1   use all available CPU cores
  --reorder     output sequences in alignment order (cleaner downstream)

After alignment this script also trims any all-gap columns and reports
basic alignment statistics (length, average pairwise identity).
"""

import os
import subprocess
import sys

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── configuration ────────────────────────────────────────────────────────────
IN_FASTA  = os.path.join("data", "ev_d68_vp1_raw.fasta")
OUT_FASTA = os.path.join("data", "ev_d68_vp1_aligned.fasta")

# Override this if mafft is not on your PATH (Windows only).
# Example: MAFFT_EXE = r"C:\Program Files\MAFFT\mafft.bat"
MAFFT_EXE = "mafft"
# ─────────────────────────────────────────────────────────────────────────────


def check_mafft():
    """Verify MAFFT is reachable, exit with a helpful message if not."""
    try:
        result = subprocess.run(
            [MAFFT_EXE, "--version"],
            capture_output=True, text=True
        )
        # MAFFT prints version to stderr
        version = result.stderr.strip() or result.stdout.strip()
        print(f"MAFFT found: {version}")
        return True
    except FileNotFoundError:
        print(
            "ERROR: MAFFT not found.\n"
            "Install it from https://mafft.cbrc.jp/alignment/software/\n"
            "Then add it to PATH, or set MAFFT_EXE in this script."
        )
        return False


def run_mafft(in_fasta, out_fasta):
    """Run MAFFT and write the aligned FASTA to out_fasta."""
    cmd = [
        MAFFT_EXE,
        "--auto",
        "--thread", "-1",
        "--reorder",
        in_fasta,
    ]
    print(f"Running: {' '.join(cmd)}")
    with open(out_fasta, "w") as fh:
        result = subprocess.run(
            cmd,
            stdout=fh,
            stderr=subprocess.PIPE,
            text=True,
        )
    if result.returncode != 0:
        print("MAFFT stderr:\n" + result.stderr)
        sys.exit(f"MAFFT exited with code {result.returncode}")

    # MAFFT progress is on stderr — show last few lines so the user knows.
    stderr_lines = result.stderr.strip().splitlines()
    if stderr_lines:
        print("MAFFT log (last 5 lines):")
        for line in stderr_lines[-5:]:
            print(f"  {line}")


def trim_all_gap_columns(aligned_fasta, out_fasta):
    """
    Remove columns where every sequence is a gap '-'.
    These arise when MAFFT inserts gaps relative to a minority of sequences
    and can artificially inflate alignment length.
    Returns the trimmed alignment.
    """
    alignment = AlignIO.read(aligned_fasta, "fasta")
    n_seq    = len(alignment)
    aln_len  = alignment.get_alignment_length()

    keep_cols = []
    for col in range(aln_len):
        column = alignment[:, col]
        if column.count("-") < n_seq:   # at least one non-gap
            keep_cols.append(col)

    trimmed_seqs = []
    for record in alignment:
        new_seq = "".join(record.seq[c] for c in keep_cols)
        trimmed_seqs.append(
            SeqRecord(Seq(new_seq), id=record.id, description="")
        )

    with open(out_fasta, "w") as fh:
        SeqIO.write(trimmed_seqs, fh, "fasta")

    removed = aln_len - len(keep_cols)
    print(f"Trimmed {removed} all-gap columns "
          f"({aln_len} → {len(keep_cols)} bp)")

    return AlignIO.read(out_fasta, "fasta")


def alignment_stats(alignment):
    """Print basic statistics about the multiple sequence alignment."""
    n_seq   = len(alignment)
    aln_len = alignment.get_alignment_length()

    # Average pairwise identity — sample at most 500 pairs for speed.
    import itertools, random
    pairs    = list(itertools.combinations(range(n_seq), 2))
    if len(pairs) > 500:
        random.seed(42)
        pairs = random.sample(pairs, 500)

    identities = []
    for i, j in pairs:
        seq_i = str(alignment[i].seq)
        seq_j = str(alignment[j].seq)
        matches = sum(a == b and a != "-" for a, b in zip(seq_i, seq_j))
        cols    = sum(a != "-" or b != "-" for a, b in zip(seq_i, seq_j))
        if cols > 0:
            identities.append(matches / cols)

    avg_id = sum(identities) / len(identities) if identities else 0.0

    print(f"\nAlignment summary:")
    print(f"  Sequences : {n_seq}")
    print(f"  Length    : {aln_len} bp")
    print(f"  Avg pairwise identity: {avg_id * 100:.1f}%")


def main():
    if not os.path.exists(IN_FASTA):
        sys.exit(f"Input not found: {IN_FASTA}\nRun fetch_sequences.py first.")

    n_seqs = sum(1 for line in open(IN_FASTA) if line.startswith(">"))
    print(f"Input: {n_seqs} sequences in {IN_FASTA}")

    if not check_mafft():
        sys.exit(1)

    tmp_aligned = OUT_FASTA + ".tmp"
    run_mafft(IN_FASTA, tmp_aligned)

    # Trim all-gap columns and overwrite the output file.
    alignment = trim_all_gap_columns(tmp_aligned, OUT_FASTA)
    os.remove(tmp_aligned)

    alignment_stats(alignment)
    print(f"\nWrote {OUT_FASTA}")


if __name__ == "__main__":
    main()
