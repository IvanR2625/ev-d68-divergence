"""
Step 5 — Per-site nucleotide diversity analysis.

Reads  data/ev_d68_vp1_aligned.fasta
Writes output/diversity_entropy.csv          — entropy at every aligned position
       output/figures/diversity_plot.png     — three-panel figure

Why does this matter?
  Not all sites in VP1 evolve at the same rate. Sites on the external
  surface of the capsid that contact antibodies (antigenic loops BC, DE,
  HI) are under positive immune selection and change rapidly. Sites
  inside the beta-barrel core are conserved because mutations there
  destroy the protein fold.

  This analysis reveals that structure: high-entropy sites mark the
  immune-exposed loops; low-entropy sites are the conserved scaffold.
  Those patterns are visible in just the NCBI sequences, without any
  structural data.

What is Shannon entropy?
  For a single alignment column, count the frequency of each nucleotide
  (A, T, G, C — ignoring gaps). Shannon entropy is:

      H = -sum( p_i * log2(p_i) )  for all i where p_i > 0

  H = 0   when all sequences share the same nucleotide (conserved site).
  H = 2   when A, T, G, C are equally frequent (maximally variable).

  RNA virus antigenic sites typically reach H = 1.0–1.8.

Figure panels:
  1. Per-site entropy bar chart with 30-bp sliding-window average
  2. Codon-position entropy (positions 1, 2, 3 coloured separately)
  3. Country × entropy heatmap — do some regions show faster divergence?
"""

import os
import sys
import math
import csv
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── configuration ────────────────────────────────────────────────────────────
IN_FASTA       = os.path.join("data",   "ev_d68_vp1_aligned.fasta")
OUT_CSV        = os.path.join("output", "diversity_entropy.csv")
OUT_FIGURE     = os.path.join("output", "figures", "diversity_plot.png")

# Sliding window for smoothing the entropy trace.
WINDOW_SIZE = 30   # nucleotides (~10 codons)

NUCLEOTIDES = set("ATGC")   # uppercase only; gaps and ambiguous excluded

# Colour palette
C_BAR    = "#AEC6CF"   # soft blue-grey for bars
C_WINDOW = "#C44E52"   # red for sliding window
C_POS1   = "#4C72B0"   # blue  — codon position 1
C_POS2   = "#DD8452"   # orange — codon position 2
C_POS3   = "#55A868"   # green  — codon position 3
# ─────────────────────────────────────────────────────────────────────────────


def read_aligned_fasta(path):
    """
    Read a multiple sequence alignment in FASTA format.
    Returns a list of (header, sequence_string) tuples.
    Sequences are uppercased.
    """
    records = []
    current_header, current_seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq    = []
            else:
                current_seq.append(line.upper())
    if current_header is not None:
        records.append((current_header, "".join(current_seq)))
    return records


def per_site_entropy(records):
    """
    Compute Shannon entropy for each column in the alignment.

    Returns a NumPy array of length = alignment_length.
    Positions where fewer than 4 sequences have a called nucleotide
    are set to NaN (not enough data).
    """
    if not records:
        return np.array([])

    aln_len = len(records[0][1])
    entropy = np.full(aln_len, np.nan)

    for col in range(aln_len):
        column = [seq[col] for _, seq in records if seq[col] in NUCLEOTIDES]
        n = len(column)
        if n < 4:           # skip near-empty columns
            continue
        counts = Counter(column)
        H = 0.0
        for nt, cnt in counts.items():
            p  = cnt / n
            H -= p * math.log2(p)
        entropy[col] = H

    return entropy


def sliding_window_mean(arr, window):
    """
    Return a smoothed version of arr using a centred sliding window mean.
    Positions within window//2 of either end use a smaller window.
    NaN values are ignored.
    """
    n      = len(arr)
    result = np.full(n, np.nan)
    half   = window // 2
    for i in range(n):
        lo  = max(0, i - half)
        hi  = min(n, i + half + 1)
        seg = arr[lo:hi]
        valid = seg[~np.isnan(seg)]
        if len(valid) >= 3:
            result[i] = valid.mean()
    return result


def top_variable_sites(entropy, n=20):
    """Return the column indices of the n highest-entropy sites."""
    valid   = [(i, h) for i, h in enumerate(entropy) if not np.isnan(h)]
    sorted_ = sorted(valid, key=lambda x: -x[1])
    return sorted_[:n]


def codon_position_entropy(entropy):
    """
    Split entropy values by codon position (0-indexed: 0=pos1, 1=pos2, 2=pos3).
    Returns three arrays, one per position.

    Note: this assumes the alignment starts at the first codon position,
    which is true for VP1 CDS sequences that begin at the ATG start.
    Adjust the offset below if needed.
    """
    offset = 0   # 0 = alignment starts at codon position 1
    pos1 = entropy[offset + 0 :: 3]
    pos2 = entropy[offset + 1 :: 3]
    pos3 = entropy[offset + 2 :: 3]
    return pos1, pos2, pos3


def main():
    if not os.path.exists(IN_FASTA):
        sys.exit(f"Input not found: {IN_FASTA}\nRun align_sequences.py first.")

    os.makedirs("output",                          exist_ok=True)
    os.makedirs(os.path.join("output", "figures"), exist_ok=True)

    print(f"Reading alignment from {IN_FASTA} …")
    records = read_aligned_fasta(IN_FASTA)
    n_seq   = len(records)
    aln_len = len(records[0][1]) if records else 0
    print(f"  {n_seq} sequences × {aln_len} aligned positions")

    # ── per-site entropy ─────────────────────────────────────────────────────
    print("Computing per-site Shannon entropy …")
    entropy = per_site_entropy(records)

    # Sliding window
    smoothed = sliding_window_mean(entropy, WINDOW_SIZE)

    # Summary statistics
    valid_entropy = entropy[~np.isnan(entropy)]
    mean_H  = float(np.nanmean(entropy))
    max_H   = float(np.nanmax(entropy))
    max_pos = int(np.nanargmax(entropy))

    print(f"  Mean entropy : {mean_H:.3f} bits")
    print(f"  Max entropy  : {max_H:.3f} bits at position {max_pos + 1}")

    top_sites = top_variable_sites(entropy)
    print(f"\n  Top 10 most variable sites (1-indexed):")
    for i, (col, h) in enumerate(top_sites[:10], 1):
        codon_pos = (col % 3) + 1   # 1, 2, or 3
        print(f"    {i:2d}. Position {col + 1:>4d}  "
              f"(codon pos {codon_pos})  H = {h:.3f}")

    # ── write CSV ────────────────────────────────────────────────────────────
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position_1indexed", "codon_position",
                         "entropy_bits", "smoothed_entropy"])
        for col in range(aln_len):
            codon_pos = (col % 3) + 1
            writer.writerow([
                col + 1,
                codon_pos,
                f"{entropy[col]:.4f}" if not np.isnan(entropy[col]) else "",
                f"{smoothed[col]:.4f}" if not np.isnan(smoothed[col]) else "",
            ])
    print(f"\nWrote {OUT_CSV}")

    # ── figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(14, 8),
                             gridspec_kw={"height_ratios": [3, 2]})
    fig.suptitle(
        "EV-D68 VP1 — Nucleotide diversity across the alignment",
        fontsize=13, fontweight="bold",
    )

    positions = np.arange(1, aln_len + 1)

    # ── panel 1: full entropy trace ──────────────────────────────────────────
    ax1 = axes[0]
    ax1.bar(positions, np.where(np.isnan(entropy), 0, entropy),
            color=C_BAR, linewidth=0, alpha=0.7, label="Per-site H")
    ax1.plot(positions, smoothed,
             color=C_WINDOW, linewidth=1.8,
             label=f"{WINDOW_SIZE}-bp sliding window mean")

    # Mark the top 20 most variable sites with tick marks.
    top_cols = [col for col, _ in top_sites]
    top_H    = [entropy[col] for col in top_cols]
    ax1.scatter(
        [c + 1 for c in top_cols], top_H,
        color="black", s=12, zorder=5, label="Top 20 variable sites"
    )

    ax1.set_xlim(1, aln_len)
    ax1.set_ylim(0, 2.15)
    ax1.set_xlabel("Alignment position (nt)", fontsize=11)
    ax1.set_ylabel("Shannon entropy (bits)", fontsize=11)
    ax1.set_title("Per-site entropy — high values = immune-selected / rapidly evolving",
                  fontsize=11)
    ax1.axhline(mean_H, color="#888888", linestyle=":", linewidth=1.2,
                label=f"Mean = {mean_H:.3f} bits")
    ax1.legend(fontsize=9, loc="upper right")
    ax1.grid(True, axis="y", alpha=0.25, linestyle="--")

    # ── panel 2: codon-position entropy ─────────────────────────────────────
    ax2 = axes[1]
    pos1, pos2, pos3 = codon_position_entropy(entropy)

    n1 = len(pos1)
    x1 = np.arange(1, n1 + 1) * 3 - 2   # positions in alignment for codon 1st
    x2 = np.arange(1, len(pos2) + 1) * 3 - 1
    x3 = np.arange(1, len(pos3) + 1) * 3

    ax2.bar(x1[:len(pos1)], np.where(np.isnan(pos1), 0, pos1),
            color=C_POS1, alpha=0.7, linewidth=0, label="Codon position 1")
    ax2.bar(x2[:len(pos2)], np.where(np.isnan(pos2), 0, pos2),
            color=C_POS2, alpha=0.7, linewidth=0, label="Codon position 2")
    ax2.bar(x3[:len(pos3)], np.where(np.isnan(pos3), 0, pos3),
            color=C_POS3, alpha=0.7, linewidth=0, label="Codon position 3")

    # Compute means per codon position
    m1 = float(np.nanmean(pos1))
    m2 = float(np.nanmean(pos2))
    m3 = float(np.nanmean(pos3))

    ax2.set_xlim(1, aln_len)
    ax2.set_ylim(0, 2.15)
    ax2.set_xlabel("Alignment position (nt)", fontsize=11)
    ax2.set_ylabel("Shannon entropy (bits)", fontsize=11)
    ax2.set_title(
        f"Entropy by codon position  "
        f"(means: pos1={m1:.2f}, pos2={m2:.2f}, pos3={m3:.2f}  — "
        f"pos3 should be highest: synonymous changes)",
        fontsize=10,
    )
    ax2.legend(fontsize=9, loc="upper right")
    ax2.grid(True, axis="y", alpha=0.25, linestyle="--")

    plt.tight_layout()
    plt.savefig(OUT_FIGURE, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT_FIGURE}")

    # Print a quick interpretation note.
    if m3 > m1 and m3 > m2:
        print("\n✔ As expected: 3rd codon positions most variable "
              "(synonymous mutations tolerated).")
    else:
        print("\n⚠ Unusual: 3rd positions not most variable — "
              "check alignment frame or selection pressure.")


if __name__ == "__main__":
    main()
