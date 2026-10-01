"""
Step 6 — Draw a coloured phylogenetic tree.

Reads  output/ev_d68_vp1.nwk   (Newick tree from build_tree.py)
       data/metadata.csv        (collection dates and countries)
Writes output/figures/phylogenetic_tree.png

Why visualise the tree?
  The ASCII tree from build_tree.py is hard to read at scale. This
  script draws a proper horizontal cladogram with:
    - Tips coloured by collection year (gradient: blue 2014 → red 2024)
    - A colour-bar legend mapping colour to year
    - Short, readable labels (year_accession instead of full header)

  The tree structure reveals whether sequences from different years
  cluster together (ancestral clades resampled) or whether there is
  clear temporal progression — more recent sequences appearing on longer
  terminal branches from older nodes. The latter pattern is expected for
  a steadily evolving virus.

Output note:
  With 150+ tips, the tree will be tall. The figure height is set
  proportional to the number of tips so every label is legible.
  For very large inputs (>300 tips) you may want to reduce LABEL_FONTSIZE.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch

from Bio import Phylo

# ── add scripts/ to path so we can import utils ──────────────────────────────
import sys; sys.path.insert(0, os.path.dirname(__file__))
from utils import load_metadata, accession_from_tip

# ── configuration ────────────────────────────────────────────────────────────
IN_NWK      = os.path.join("output", "ev_d68_vp1.nwk")
IN_METADATA = os.path.join("data",   "metadata.csv")
OUT_FIGURE  = os.path.join("output", "figures", "phylogenetic_tree.png")

LABEL_FONTSIZE = 6      # reduce if you have many sequences
FIGURE_WIDTH   = 12     # inches
PTS_PER_TIP    = 12     # vertical pixels-equivalent per tip (controls height)
COLORMAP       = "coolwarm"   # blue (early) → red (recent)
# ─────────────────────────────────────────────────────────────────────────────


def rename_tips(tree, meta):
    """
    Replace full FASTA headers with short labels: 'YYYY_ACCCODE'.
    Labels longer than this clutter the tree.
    Returns a dict { new_name: year } for colouring.
    """
    name_to_year = {}
    for tip in tree.get_terminals():
        accession = accession_from_tip(tip.name)
        year      = meta.get(accession, {}).get("year", None)
        short_acc = accession[:10]   # first 10 chars of accession
        new_name  = f"{year}_{short_acc}" if year else short_acc
        tip.name  = new_name
        name_to_year[new_name] = year
    return name_to_year


def year_color_map(name_to_year, cmap_name=COLORMAP):
    """
    Return a dict { tip_name: hex_color } mapping each tip to a colour
    from the chosen colormap, scaled over the observed year range.
    """
    years  = [y for y in name_to_year.values() if y is not None]
    if not years:
        return {k: "#888888" for k in name_to_year}

    y_min, y_max = min(years), max(years)
    cmap  = cm.get_cmap(cmap_name)
    norm  = mcolors.Normalize(vmin=y_min, vmax=y_max)

    color_map = {}
    for name, year in name_to_year.items():
        if year is None:
            color_map[name] = "#888888"
        else:
            rgba = cmap(norm(year))
            color_map[name] = mcolors.to_hex(rgba)
    return color_map, norm, cmap, y_min, y_max


def draw_tree(tree, color_map, norm, cmap, y_min, y_max):
    """
    Draw the phylogenetic tree using Biopython's Phylo.draw().
    Tips are coloured by year.
    """
    n_tips  = len(tree.get_terminals())
    fig_h   = max(8, n_tips * PTS_PER_TIP / 72)   # 72 pt/inch
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, fig_h))

    # Biopython draw() colours tips by passing a label_colors dict.
    label_colors = color_map   # { tip_name: hex_color }

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_colors=label_colors,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )

    # Style axes
    ax.set_title(
        f"EV-D68 VP1 neighbour-joining tree  ({n_tips} sequences)\n"
        "Tips coloured by collection year  (blue = earlier, red = later)",
        fontsize=11, fontweight="bold",
    )
    ax.set_xlabel("Substitutions per site", fontsize=10)

    # Fix tiny default font on tip labels
    for text_obj in ax.texts:
        text_obj.set_fontsize(LABEL_FONTSIZE)

    # Colorbar legend
    sm  = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical",
                        fraction=0.015, pad=0.01)
    cbar.set_label("Collection year", fontsize=10)
    tick_years = list(range(y_min, y_max + 1, 2))
    cbar.set_ticks(tick_years)
    cbar.ax.set_yticklabels([str(y) for y in tick_years], fontsize=9)

    return fig


def main():
    for path in (IN_NWK, IN_METADATA):
        if not os.path.exists(path):
            sys.exit(f"Input not found: {path}\nRun build_tree.py first.")

    os.makedirs(os.path.join("output", "figures"), exist_ok=True)

    print(f"Reading tree from {IN_NWK} …")
    tree = Phylo.read(IN_NWK, "newick")
    n_tips = len(tree.get_terminals())
    print(f"  {n_tips} tips")

    print(f"Reading metadata from {IN_METADATA} …")
    meta = load_metadata(IN_METADATA)

    # Rename tips to short labels and build colour mapping
    name_to_year = rename_tips(tree, meta)
    color_map, norm, cmap, y_min, y_max = year_color_map(name_to_year)

    # Count tips that got a year colour vs. grey
    coloured = sum(1 for y in name_to_year.values() if y is not None)
    print(f"  {coloured} tips coloured by year, "
          f"{n_tips - coloured} without date (grey)")

    print("Drawing tree …")
    fig = draw_tree(tree, color_map, norm, cmap, y_min, y_max)

    plt.tight_layout()
    plt.savefig(OUT_FIGURE, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {OUT_FIGURE}")

    # Print a brief interpretation hint
    print(
        "\nWhat to look for in the tree:"
        "\n  - Tips of the same colour (year) should cluster if sampling"
        " is sparse, but with dense sampling they may spread throughout the tree."
        "\n  - Long terminal branches on recent (red) tips = fast-evolving lineages."
        "\n  - A 'broom' pattern (many short tips from recent years) suggests"
        " a population expansion."
    )


if __name__ == "__main__":
    main()
