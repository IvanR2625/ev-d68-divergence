"""
Step 4 — Root-to-tip divergence plot with clock-outlier detection.

Reads  output/ev_d68_vp1.nwk      (midpoint-rooted Newick tree)
       data/metadata.csv           (accession, date, country)
Writes output/figures/divergence_clock.png
       output/outliers.csv         (sequences deviating > OUTLIER_SD from
                                    the regression line)

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

What is a clock outlier?
  A sequence is flagged as an outlier when its root-to-tip distance
  deviates more than OUTLIER_SD standard deviations from the regression
  line (i.e., its residual is unusually large or small).

  Possible explanations:
    - The collection date in the GenBank record is wrong (misdated sample).
    - The sequence is from a divergent lineage or contains extra mutations.
    - Lab contamination or sequencing artefact.

  Outliers do NOT invalidate the clock — they reveal sequences worth
  inspecting. After running this step, open output/outliers.csv and
  look up each accession on NCBI to see whether the date looks plausible.

Figure panels:
  1. Root-to-tip scatter + regression; outliers highlighted in orange.
  2. Sequences per year (bar chart).
  3. Divergence distribution per year (box plot).
"""

import os
import sys
import csv
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats

from Bio import Phylo

# ── add scripts/ to path so we can import utils ──────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from utils import load_metadata, decimal_year, dec_to_dt, accession_from_tip

# ── configuration ────────────────────────────────────────────────────────────
IN_NWK       = os.path.join("output", "ev_d68_vp1.nwk")
IN_METADATA  = os.path.join("data",   "metadata.csv")
OUT_FIGURE   = os.path.join("output", "figures", "divergence_clock.png")
OUT_OUTLIERS = os.path.join("output", "outliers.csv")

VIRUS_NAME   = "Enterovirus D68"
GENE_NAME    = "VP1"

# Residuals beyond this many SD from the regression are flagged as outliers.
OUTLIER_SD   = 2.0

PALETTE = {
    "scatter":  "#4C72B0",   # muted blue — normal sequences
    "outlier":  "#FF7C00",   # orange — outliers
    "regline":  "#C44E52",   # muted red — regression line
    "ci":       "#C44E52",   # same but transparent for confidence band
    "bar":      "#55A868",   # muted green
    "violin":   "#8172B2",   # muted purple
}
# ─────────────────────────────────────────────────────────────────────────────


def root_to_tip_distances(tree):
    """Return { terminal_name: root_to_tip_distance } for a rooted tree."""
    root    = tree.root
    tip_map = {}
    for tip in tree.get_terminals():
        tip_map[tip.name] = tree.distance(root, tip)
    return tip_map


def detect_outliers(x, y, slope, intercept, n_sd=OUTLIER_SD):
    """
    Return a boolean mask: True where |residual| > n_sd * std(residuals).
    """
    residuals  = y - (slope * x + intercept)
    std_resid  = np.std(residuals)
    return np.abs(residuals) > n_sd * std_resid, residuals


def regression_ci(x, y, slope, intercept, confidence=0.95):
    """
    Compute pointwise confidence band for the regression line.

    Returns (x_plot, y_lo, y_hi) where y_lo / y_hi are the lower and
    upper bounds of the confidence interval at each x_plot value.
    Uses the standard formula for a simple linear regression CI.
    """
    n     = len(x)
    x_bar = x.mean()
    se_sq = np.sum((y - (slope * x + intercept))**2) / (n - 2)
    x_plt = np.linspace(x.min(), x.max(), 200)
    # Variance of the fitted value at each x_plt point
    var_fit = se_sq * (1/n + (x_plt - x_bar)**2 / np.sum((x - x_bar)**2))

    from scipy.stats import t as t_dist
    t_crit  = t_dist.ppf((1 + confidence) / 2, df=n - 2)
    margin  = t_crit * np.sqrt(var_fit)

    y_hat = slope * x_plt + intercept
    return x_plt, y_hat - margin, y_hat + margin


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

    # ── root-to-tip distances ────────────────────────────────────────────────
    rtd_map = root_to_tip_distances(tree)

    dates_dec  = []
    dates_dt   = []
    divs       = []
    years      = []
    accessions = []

    skipped = 0
    for tip_name, rtd in rtd_map.items():
        acc = accession_from_tip(tip_name)
        if acc not in meta:
            skipped += 1
            continue
        m = meta[acc]
        dates_dec.append(decimal_year(m["date"]))
        dates_dt.append(m["date"])
        divs.append(rtd)
        years.append(m["year"])
        accessions.append(acc)

    if skipped:
        print(f"  {skipped} tips had no metadata match.")
    print(f"  {len(dates_dec)} tip–date pairs for plotting.")

    if len(dates_dec) < 3:
        sys.exit("Too few matched sequences. Check accession IDs.")

    x = np.array(dates_dec)
    y = np.array(divs)

    # ── linear regression ────────────────────────────────────────────────────
    slope, intercept, r_value, p_value, se = stats.linregress(x, y)
    r2 = r_value ** 2
    x_line     = np.linspace(x.min(), x.max(), 200)
    y_line     = slope * x_line + intercept
    x_line_dt  = [dec_to_dt(v) for v in x_line]

    print(f"\nRegression results:")
    print(f"  Slope (evolutionary rate): {slope:.3e} subs/site/year")
    print(f"  R²  : {r2:.3f}")
    print(f"  p   : {p_value:.2e}")

    # ── outlier detection ────────────────────────────────────────────────────
    outlier_mask, residuals = detect_outliers(x, y, slope, intercept)
    n_outliers = outlier_mask.sum()
    print(f"\nClock outliers (|residual| > {OUTLIER_SD} SD): {n_outliers}")

    if n_outliers > 0:
        outlier_rows = []
        for i, is_out in enumerate(outlier_mask):
            if is_out:
                direction = "above" if residuals[i] > 0 else "below"
                print(f"  {accessions[i]}  {years[i]}  "
                      f"resid={residuals[i]:+.4f}  ({direction} line)")
                outlier_rows.append({
                    "accession":   accessions[i],
                    "year":        years[i],
                    "date":        dates_dt[i].strftime("%Y-%m-%d"),
                    "divergence":  f"{divs[i]:.5f}",
                    "residual":    f"{residuals[i]:+.5f}",
                    "sd_from_line": f"{abs(residuals[i]) / np.std(residuals):.2f}",
                    "direction":   direction,
                })
        with open(OUT_OUTLIERS, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(outlier_rows[0].keys()))
            w.writeheader()
            w.writerows(outlier_rows)
        print(f"  Saved → {OUT_OUTLIERS}")
    else:
        print("  No outliers detected — strong clock signal.")

    # ── 95% confidence band ──────────────────────────────────────────────────
    x_ci, y_lo, y_hi = regression_ci(x, y, slope, intercept)
    x_ci_dt = [dec_to_dt(v) for v in x_ci]

    # ── figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        f"{VIRUS_NAME} {GENE_NAME}  |  Longitudinal divergence  "
        f"({int(x.min())}–{int(x.max())})",
        fontsize=13, fontweight="bold", y=1.01,
    )

    # ── panel 1: clock scatter ───────────────────────────────────────────────
    ax1 = axes[0]

    # Normal sequences (not outliers)
    normal_mask = ~outlier_mask
    ax1.scatter(
        np.array(dates_dt)[normal_mask],
        y[normal_mask],
        color=PALETTE["scatter"], alpha=0.65, s=18, linewidths=0,
        label=f"Sequence (n={normal_mask.sum()})",
        zorder=3,
    )

    # Outliers
    if n_outliers:
        ax1.scatter(
            np.array(dates_dt)[outlier_mask],
            y[outlier_mask],
            color=PALETTE["outlier"], alpha=0.9, s=35,
            marker="D", linewidths=0.5, edgecolors="white",
            label=f"Clock outlier (n={n_outliers})",
            zorder=4,
        )

    # 95% confidence band
    ax1.fill_between(
        x_ci_dt, y_lo, y_hi,
        color=PALETTE["ci"], alpha=0.12, label="95% CI",
    )

    # Regression line
    ax1.plot(
        x_line_dt, y_line,
        color=PALETTE["regline"], linewidth=2,
        label=f"Rate = {slope:.2e} subs/site/yr",
        zorder=5,
    )

    ax1.set_xlabel("Collection date", fontsize=11)
    ax1.set_ylabel("Root-to-tip divergence (subs/site)", fontsize=11)
    ax1.set_title("Molecular clock", fontsize=12)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.xaxis.set_major_locator(mdates.YearLocator(2))
    plt.setp(ax1.get_xticklabels(), rotation=30, ha="right")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(True, alpha=0.25, linestyle="--")

    ax1.annotate(
        f"$R^2 = {r2:.3f}$\n$p = {p_value:.1e}$",
        xy=(0.97, 0.05), xycoords="axes fraction",
        fontsize=10, va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
    )

    # ── panel 2: sequences per year ──────────────────────────────────────────
    ax2 = axes[1]
    year_list   = sorted(set(years))
    year_counts = [years.count(yr) for yr in year_list]
    # Colour bars by outlier fraction
    outlier_by_year = {yr: 0 for yr in year_list}
    for i, is_out in enumerate(outlier_mask):
        if is_out:
            outlier_by_year[years[i]] += 1

    bars = ax2.bar(year_list, year_counts,
                   color=PALETTE["bar"], edgecolor="white", linewidth=0.5)
    # Overlay outlier counts in orange
    outlier_counts = [outlier_by_year[yr] for yr in year_list]
    ax2.bar(year_list, outlier_counts,
            color=PALETTE["outlier"], edgecolor="white", linewidth=0.5,
            alpha=0.85, label="Outliers")

    ax2.set_xlabel("Year", fontsize=11)
    ax2.set_ylabel("Number of sequences", fontsize=11)
    ax2.set_title("Sequences per year", fontsize=12)
    ax2.set_xticks(year_list)
    ax2.set_xticklabels(year_list, rotation=45, ha="right")
    ax2.grid(True, axis="y", alpha=0.25, linestyle="--")
    if any(outlier_counts):
        ax2.legend(fontsize=9)
    for bar, count in zip(bars, year_counts):
        ax2.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + max(year_counts) * 0.01,
                 str(count), ha="center", va="bottom", fontsize=8)

    # ── panel 3: divergence by year (box) ────────────────────────────────────
    ax3 = axes[2]
    data_by_year = {yr: [] for yr in year_list}
    for yr, div in zip(years, divs):
        data_by_year[yr].append(div)

    ax3.boxplot(
        [data_by_year[yr] for yr in year_list],
        labels=year_list,
        patch_artist=True,
        medianprops=dict(color=PALETTE["regline"], linewidth=2),
        boxprops=dict(facecolor=PALETTE["violin"], alpha=0.5),
        whiskerprops=dict(linestyle="--"),
        flierprops=dict(marker="o", markerfacecolor=PALETTE["outlier"],
                        markersize=4, alpha=0.6),
    )
    ax3.set_xlabel("Year", fontsize=11)
    ax3.set_ylabel("Root-to-tip divergence (subs/site)", fontsize=11)
    ax3.set_title("Divergence distribution per year", fontsize=12)
    ax3.set_xticklabels(year_list, rotation=45, ha="right")
    ax3.grid(True, axis="y", alpha=0.25, linestyle="--")

    plt.tight_layout()
    plt.savefig(OUT_FIGURE, dpi=150, bbox_inches="tight")
    print(f"\nWrote {OUT_FIGURE}")

    # ── plain-text summary (also captured by generate_report.py) ─────────────
    print("\n── Summary ──────────────────────────────────────────────────────────")
    print(f"Virus              : {VIRUS_NAME}")
    print(f"Gene               : {GENE_NAME}")
    print(f"Sequences          : {len(dates_dec)}")
    print(f"Date range         : {min(dates_dt).strftime('%Y-%m-%d')} – "
          f"{max(dates_dt).strftime('%Y-%m-%d')}")
    print(f"Evo. rate          : {slope:.3e} substitutions/site/year")
    print(f"R²                 : {r2:.3f}")
    print(f"p-value            : {p_value:.2e}")
    print(f"Mean divergence    : {np.mean(y):.4f} subs/site")
    print(f"Clock outliers     : {n_outliers} "
          f"({'see ' + OUT_OUTLIERS if n_outliers else 'none detected'})")
    print("──────────────────────────────────────────────────────────────────────")


if __name__ == "__main__":
    main()
