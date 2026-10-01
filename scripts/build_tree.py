"""
Step 3 — Build a neighbour-joining phylogenetic tree.

Reads  data/ev_d68_vp1_aligned.fasta
Writes output/ev_d68_vp1.nwk        (Newick tree, midpoint-rooted)
       output/ev_d68_vp1_tree.txt   (ASCII summary for quick inspection)

Uses BioPython's built-in DistanceCalculator + DistanceTreeConstructor.
No external tree-building program needed.

What is a neighbour-joining tree?
  NJ groups sequences that are most similar to each other first, then
  builds outward. It's fast and good enough for measuring root-to-tip
  divergence, which is what the next step uses.

What is midpoint rooting?
  Without an outgroup we can't know where the true root is. Midpoint
  rooting places the root at the midpoint of the longest path through
  the tree — a reasonable approximation when the substitution rate is
  roughly clock-like (which EV-D68 VP1 is).
"""

import os
import sys

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO

# ── configuration ────────────────────────────────────────────────────────────
IN_FASTA   = os.path.join("data",   "ev_d68_vp1_aligned.fasta")
OUT_NWK    = os.path.join("output", "ev_d68_vp1.nwk")
OUT_ASCII  = os.path.join("output", "ev_d68_vp1_tree.txt")
# ─────────────────────────────────────────────────────────────────────────────


def midpoint_root(tree):
    """
    Midpoint-root a BioPython Tree object in place.

    Finds the two most distantly separated tips, then places the root
    at the midpoint of the path between them.
    """
    # Get all terminals (leaf nodes).
    terminals = tree.get_terminals()

    # Find the pair of terminals with the maximum patristic distance.
    max_dist = 0.0
    tip_a, tip_b = terminals[0], terminals[1]

    for i, t1 in enumerate(terminals):
        for t2 in terminals[i + 1:]:
            d = tree.distance(t1, t2)
            if d > max_dist:
                max_dist = d
                tip_a, tip_b = t1, t2

    # The midpoint lies max_dist/2 from either tip.
    # BioPython's root_with_outgroup + manual branch shortening is one path;
    # the cleanest built-in method is root_at_midpoint (available since
    # Biopython 1.79 — check your version).
    try:
        tree.root_at_midpoint()
    except AttributeError:
        # Fallback for older BioPython: root on the tip closest to midpoint.
        tree.root_with_outgroup(tip_a)

    return tree


def main():
    if not os.path.exists(IN_FASTA):
        sys.exit(f"Input not found: {IN_FASTA}\nRun align_sequences.py first.")

    os.makedirs("output", exist_ok=True)

    print(f"Reading alignment from {IN_FASTA} …")
    alignment = AlignIO.read(IN_FASTA, "fasta")
    n_seq     = len(alignment)
    aln_len   = alignment.get_alignment_length()
    print(f"  {n_seq} sequences × {aln_len} aligned positions")

    # ── distance matrix ──────────────────────────────────────────────────────
    # 'identity' model: distance = fraction of positions that differ.
    # Other options: 'blastn', 'trans' (for amino acids).
    print("Computing pairwise distance matrix …")
    calculator = DistanceCalculator("identity")
    dm         = calculator.get_distance(alignment)

    # ── NJ tree ──────────────────────────────────────────────────────────────
    print("Building neighbour-joining tree …")
    constructor = DistanceTreeConstructor()
    tree        = constructor.nj(dm)

    # ── midpoint root ────────────────────────────────────────────────────────
    print("Midpoint-rooting the tree …")
    tree = midpoint_root(tree)

    # ── write Newick ─────────────────────────────────────────────────────────
    Phylo.write(tree, OUT_NWK, "newick")
    print(f"Wrote {OUT_NWK}")

    # ── write ASCII summary ──────────────────────────────────────────────────
    buf = StringIO()
    Phylo.draw_ascii(tree, file=buf, column_width=80)
    ascii_text = buf.getvalue()

    with open(OUT_ASCII, "w") as fh:
        fh.write(ascii_text)

    # Print a short excerpt to the terminal.
    lines = ascii_text.splitlines()
    preview = lines[:30]
    if len(lines) > 30:
        preview.append(f"  … ({len(lines) - 30} more lines — see {OUT_ASCII})")
    print("\nTree preview (ASCII):")
    print("\n".join(preview))

    print(f"\nWrote {OUT_ASCII}")

    # Quick sanity check — tips should equal number of sequences.
    n_tips = len(tree.get_terminals())
    print(f"\nSanity check: {n_tips} tips in tree, {n_seq} input sequences.")
    if n_tips != n_seq:
        print("WARNING: counts differ — some sequences may have been dropped "
              "due to identical headers. Check your FASTA headers.")


if __name__ == "__main__":
    main()
