"""
Step 4 — Root-to-tip divergence plot.

Reads  output/ev_d68_vp1.nwk      (midpoint-rooted Newick tree)
       data/metadata.csv           (accession, date, country)
Writes output/figures/divergence_clock.png

What is root-to-tip divergence?
  For each leaf (sampled sequence), we measure the total branch length
  from the root to that leaf — i.e., the number of substitutions per
  site accumulated since the common ancestor. If the virus evolves at
  a roughly constant rate (a "molecular clock"), this number should
  increase linearly with time. A scatter plot of divergence vs.
  collection date and a regression line tests that expectation.

  The slope of the regression is the evolutionary rate:
    slope ≈ substitutions per site per year

  The R² value tells you how clock-like the data are (R² > 0.8 is
  typical for RNA viruses with good sampling).

Outputs three panels in a single figure:
  1. Root-to-tip divergence vs. collection date (scatter + regression)
  2. Sequences per year (bar chart)
  3. Divergence distribution by year (violin/box)
"""

import os
import sys
import csv
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")          # non-interactive backend — safe on any system
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats

from Bio import Phylo

# ── configuration ────────────────────────────────────────────────────────────
IN_NWK      = os.path.join("output", "ev_d68_vp1.nwk")
IN_METADATA = os.path.join("data",   "metadata.csv")
OUT_FIGURE  = os.path.join("output", "figures", "divergence_clock.png")

VIRUS_NAME  = "Enterovirus D68"
GENE_NAME   = "VP1"

# Colour palette — feels scientific, prints in greyscale reasonably well.
PALETTE = {
    "scatter": "#4C72B0",   # muted blue
    "regline": "#C44E52",   # muted red
    "bar":     "#55A868",   # muted green
    "violin":  "#8172B2",   # muted purple
}
# ─────────────────────────────────────────────────────────────────────────────


def root_to_tip_distances(tree):
    """
    Return a dict  { terminal_name: root_to_tip_distance }
    by summing branch lengths from the root to every leaf.

    BioPython's tree.distance(a, b) gives the path between any two
    nodes; here we use the root as node a.
    """
    root    = tree.root
    tip_map = {}
    for tip in tree.get_terminals():
        tip_map[tip.name] = tree.distance(root, tip)
    return tip_map


def load_metadata(csv_path):
    """
    Return a dict  { accession: {'date': datetime, 'year': int, ...} }
    """
    meta = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                date_obj = datetime.strptime(row["date"], "%Y-%m-%d")
            except ValueError:
                continue
            meta[row["accession"]] = {
                "date":    date_obj,
                "year":    int(row["year"]),
                "country": row["country"],
            }
    return meta


def decimal_year(dt):
    """Convert a datetime to a decimal year (e.g. 2014-07-01 → 2014.5)."""
    year_start = datetime(dt.year, 1, 1)
    year_end   = datetime(dt.year + 1, 1, 1)
    fraction   = (dt - year_start) / (year_end - year_start)
    return dt.year + fraction


def parse_accession_from_tip(tip_name):
    """
    Tip names in the Newick tree come from FASTA headers:
      accession|date|country
    Return just the accession part.
    """
    return tip_name.split("|")[0]


def main():
    for path in (IN_NWK, IN_METADATA):
        if not os.path.exists(path):
            sys.exit(f"Input not found: {path}\nRun the previous steps first.")

    os.makedirs(os.path.join("output", "figures"), exist_ok=True)

    # ── load tree and metadata ───────────────────────────────────────────────
    print(f"Reading tree from {IN_NWK} …")
    tree = Phylo.read(IN_NWK, "newick")

    print(f"Reading metadata from {IN_METADATA} …")
    meta = load_metadata(IN_METADATA)

    # ── compute root-to-tip distances ────────────────────────────────────────
    rtd_map = root_to_tip_distances(tree)

    # Match tip names to metadata.
    dates      = []   # decimal year
    dates_dt   = []   # datetime objects (for x-axis ticks)
    divergences = []
    years      = []

    skipped = 0
    for tip_name, rtd in rtd_map.items():
        accession = parse_accession_from_tip(tip_name)
        if accession not in meta:
            skipped += 1
            continue
        m = meta[accession]
        dates.append(decimal_year(m["date"]))
        dates_dt.append(m["date"])
        divergences.append(rtd)
        years.append(m["year"])

    if skipped:
        print(f"  {skipped} tips had no metadata match (check accession format).")
    print(f"  {len(dates)} tip–date pairs for plotting.")

    if len(dates) < 3:
        sys.exit("Too few matched sequences to plot. Check that accession IDs "
                 "in the tree match those in metadata.csv.")

    x = np.array(dates)
    y = np.array(divergences)

    # ── linear regression ────────────────────────────────────────────────────
    slope, intercept, r_value, p_value, se = stats.linregress(x, y)
    r2  = r_value ** 2
    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = slope * x_line + intercept

    rate_per_year = slope   # substitutions per site per year
    rate_per_site_per_day = slope / 365.25

    print(f"\nRegression results:")
    print(f"  Slope (evolutionary rate): {rate_per_year:.2e} subs/site/year")
    print(f"  R²  : {r2:.3f}")
    print(f"  p   : {p_value:.2e}")

    # ── figure: three-panel layout ───────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        f"{VIRUS_NAME} {GENE_NAME}  |  Longitudinal divergence  "
        f"({int(x.min())}–{int(x.max())})",
        fontsize=13, fontweight="bold", y=1.01,
    )

    # ── panel 1: root-to-tip scatter + regression ────────────────────────────
    ax1 = axes[0]
    ax1.scatter(dates_dt, y,
                color=PALETTE["scatter"], alpha=0.6, s=18, linewidths=0,
                label="Sequence")

    # Convert decimal year back to datetime for the regression line.
    def dec_to_dt(dy):
        yr    = int(dy)
        frac  = dy - yr
        days  = (datetime(yr + 1, 1, 1) - datetime(yr, 1, 1)).days
        return datetime(yr, 1, 1) + __import__("datetime").timedelta(days=frac * days)

    x_line_dt = [dec_to_dt(v) for v in x_line]
    ax1.plot(x_line_dt, y_line,
             color=PALETTE["regline"], linewidth=2,
             label=f"Rate = {rate_per_year:.2e} subs/site/yr\n$R^2$ = {r2:.3f}")

    ax1.set_xlabel("Collection date", fontsize=11)
    ax1.set_ylabel("Root-to-tip divergence (subs/site)", fontsize=11)
    ax1.set_title("Molecular clock", fontsize=12)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.xaxis.set_major_locator(mdates.YearLocator(2))
    plt.setp(ax1.get_xticklabels(), rotation=30, ha="right")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, linestyle="--")

    # Annotate R² prominently.
    ax1.annotate(
        f"$R^2 = {r2:.2f}$\n$p = {p_value:.1e}$",
        xy=(0.05, 0.93), xycoords="axes fraction",
        fontsize=10, va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )

    # ── panel 2: sequences per year (bar) ────────────────────────────────────
    ax2 = axes[1]
    year_list  = sorted(set(years))
    year_counts = [years.count(yr) for yr in year_list]
    ax2.bar(year_list, year_counts,
            color=PALETTE["bar"], edgecolor="white", linewidth=0.5)
    ax2.set_xlabel("Year", fontsize=11)
    ax2.set_ylabel("Number of sequences", fontsize=11)
    ax2.set_title("Sequences per year", fontsize=12)
    ax2.set_xticks(year_list)
    ax2.set_xticklabels(year_list, rotation=45, ha="right")
    ax2.grid(True, axis="y", alpha=0.3, linestyle="--")
    for bar, count in zip(ax2.patches, year_counts):
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 0.3,
                 str(count), ha="center", va="bottom", fontsize=8)

    # ── panel 3: divergence by year (box plot) ───────────────────────────────
    ax3 = axes[2]
    data_by_year = {yr: [] for yr in year_list}
    for yr, div in zip(years, divergences):
        data_by_year[yr].append(div)

    bp = ax3.boxplot(
        [data_by_year[yr] for yr in year_list],
        labels=year_list,
        patch_artist=True,
        medianprops=dict(color=PALETTE["regline"], linewidth=2),
        boxprops=dict(facecolor=PALETTE["violin"], alpha=0.5),
        whiskerprops=dict(linestyle="--"),
    )
    ax3.set_xlabel("Year", fontsize=11)
    ax3.set_ylabel("Root-to-tip divergence (subs/site)", fontsize=11)
    ax3.set_title("Divergence distribution per year", fontsize=12)
    ax3.set_xticklabels(year_list, rotation=45, ha="right")
    ax3.grid(True, axis="y", alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(OUT_FIGURE, dpi=150, bbox_inches="tight")
    print(f"\nWrote {OUT_FIGURE}")

    # ── print a plain-text summary for the README / methods section ──────────
    print("\n── Summary ─────────────────────────────────────────────────────────")
    print(f"Virus          : {VIRUS_NAME}")
    print(f"Gene           : {GENE_NAME}")
    print(f"Sequences      : {len(dates)}")
    print(f"Date range     : {min(dates_dt).strftime('%Y-%m-%d')} – "
          f"{max(dates_dt).strftime('%Y-%m-%d')}")
    print(f"Evo. rate      : {rate_per_year:.3e} substitutions/site/year")
    print(f"R²             : {r2:.3f}")
    print(f"p-value        : {p_value:.2e}")
    print(f"Mean divergence: {np.mean(y):.4f} subs/site")
    print("────────────────────────────────────────────────────────────────────")


if __name__ == "__main__":
    main()
