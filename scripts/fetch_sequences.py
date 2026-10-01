"""
Step 7 — Generate a Markdown summary report.

Reads  output/outliers.csv             (if it exists)
       output/diversity_entropy.csv    (if it exists)
       data/metadata.csv
       output/ev_d68_vp1.nwk
Writes output/REPORT.md

The report collects all key numbers from the pipeline in one place,
links to the figures, explains what each result means, and calls out
anything worth investigating (outliers, unexpected codon-position
entropy, low R², etc.).

It is designed to paste directly into a science fair write-up or GitHub
README, and to be a record you can hand to a mentor.
"""

import os
import sys
import csv
import math
from datetime import datetime

import numpy as np
from scipy import stats
from Bio import Phylo, AlignIO

# ── add scripts/ to path so we can import utils ──────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from utils import load_metadata, decimal_year, accession_from_tip

# ── configuration ────────────────────────────────────────────────────────────
IN_NWK       = os.path.join("output", "ev_d68_vp1.nwk")
IN_METADATA  = os.path.join("data",   "metadata.csv")
IN_OUTLIERS  = os.path.join("output", "outliers.csv")
IN_ENTROPY   = os.path.join("output", "diversity_entropy.csv")
IN_ALIGNED   = os.path.join("data",   "ev_d68_vp1_aligned.fasta")
OUT_REPORT   = os.path.join("output", "REPORT.md")

VIRUS_NAME   = "Enterovirus D68"
GENE_NAME    = "VP1"
OUTLIER_SD   = 2.0
# ─────────────────────────────────────────────────────────────────────────────


def load_outliers(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_entropy_csv(path):
    """Return (positions, entropy values, smoothed values) as numpy arrays."""
    if not os.path.exists(path):
        return None, None, None
    pos, H, Hs = [], [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            pos.append(int(row["position_1indexed"]))
            H.append(float(row["entropy_bits"]) if row["entropy_bits"] else float("nan"))
            Hs.append(float(row["smoothed_entropy"]) if row["smoothed_entropy"] else float("nan"))
    return np.array(pos), np.array(H), np.array(Hs)


def root_to_tip_distances(tree):
    root = tree.root
    return {tip.name: tree.distance(root, tip) for tip in tree.get_terminals()}


def compute_clock_stats(tree, meta):
    """Return (slope, intercept, r2, p, x_arr, y_arr, dates_dt_arr)."""
    rtd_map = root_to_tip_distances(tree)
    x, y, dates_dt = [], [], []
    for tip_name, rtd in rtd_map.items():
        acc = accession_from_tip(tip_name)
        if acc not in meta:
            continue
        m = meta[acc]
        x.append(decimal_year(m["date"]))
        y.append(rtd)
        dates_dt.append(m["date"])
    x, y = np.array(x), np.array(y)
    slope, intercept, r_val, p, _ = stats.linregress(x, y)
    return slope, intercept, r_val**2, p, x, y, dates_dt


def interpret_r2(r2):
    if r2 >= 0.85:
        return "strong (≥ 0.85) — the virus shows clear clock-like evolution"
    elif r2 >= 0.70:
        return "moderate (0.70–0.84) — reasonable clock signal; adding more sequences per year would strengthen it"
    elif r2 >= 0.50:
        return "weak (0.50–0.69) — consider removing outliers or increasing sampling density"
    else:
        return "very weak (< 0.50) — check for sampling gaps, misdated records, or strong positive selection"


def compare_to_literature(rate):
    """Compare estimated rate to published EV-D68 VP1 estimates."""
    lit_lo, lit_hi = 3e-3, 7e-3
    if lit_lo <= rate <= lit_hi:
        status = "✅ within the published range"
    elif rate < lit_lo:
        status = f"⚠ below the published range — may reflect insufficient sampling or shorter sequences"
    else:
        status = f"⚠ above the published range — possibly inflated by recombination or hypervariable samples"
    return status, lit_lo, lit_hi


def top_entropy_sites(H_arr, n=5):
    """Return (1-indexed position, entropy) for the top n sites."""
    valid = [(i + 1, H_arr[i]) for i in range(len(H_arr)) if not math.isnan(H_arr[i])]
    return sorted(valid, key=lambda t: -t[1])[:n]


def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M")


def section(title, level=2):
    prefix = "#" * level
    return f"\n{prefix} {title}\n"


def main():
    for path in (IN_NWK, IN_METADATA):
        if not os.path.exists(path):
            sys.exit(f"Required input missing: {path}\nRun steps 1–4 first.")

    os.makedirs("output", exist_ok=True)

    print("Loading data for report …")

    meta     = load_metadata(IN_METADATA)
    tree     = Phylo.read(IN_NWK, "newick")
    outliers = load_outliers(IN_OUTLIERS)
    pos_arr, H_arr, Hs_arr = load_entropy_csv(IN_ENTROPY)

    slope, intercept, r2, p_val, x, y, dates_dt = compute_clock_stats(tree, meta)
    rate_status, lit_lo, lit_hi = compare_to_literature(slope)

    n_seq     = len(x)
    year_min  = int(x.min())
    year_max  = int(x.max())
    date_min  = min(dates_dt).strftime("%Y-%m-%d")
    date_max  = max(dates_dt).strftime("%Y-%m-%d")
    mean_div  = float(np.mean(y))
    n_tips    = len(tree.get_terminals())

    countries = sorted(set(m["country"] for m in meta.values()))
    n_countries = len(countries)

    # Alignment stats (if available)
    aln_len = None
    if os.path.exists(IN_ALIGNED):
        try:
            aln = AlignIO.read(IN_ALIGNED, "fasta")
            aln_len = aln.get_alignment_length()
        except Exception:
            pass

    # ── build the Markdown document ──────────────────────────────────────────
    lines = []

    lines.append("# EV-D68 VP1 Divergence Analysis — Results Report")
    lines.append(f"\n*Generated: {now_str()}*\n")
    lines.append(
        "> **Disclaimer:** This analysis is for educational and research "
        "purposes. Not validated for clinical use.\n"
    )

    # ── 1. Dataset ───────────────────────────────────────────────────────────
    lines.append(section("1. Dataset"))
    lines.append(f"| Property | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Virus | {VIRUS_NAME} |")
    lines.append(f"| Gene | {GENE_NAME} (VP1 capsid protein) |")
    lines.append(f"| Source | NCBI Nucleotide (GenBank) |")
    lines.append(f"| Sequences in tree | {n_tips} |")
    lines.append(f"| Sequences with dates | {n_seq} |")
    lines.append(f"| Date range | {date_min} – {date_max} ({year_min}–{year_max}) |")
    lines.append(f"| Countries represented | {n_countries} |")
    if aln_len:
        lines.append(f"| Alignment length | {aln_len} bp |")
    lines.append("")

    # Country list
    if countries:
        lines.append(f"**Countries:** {', '.join(countries[:20])}"
                     + (" …" if len(countries) > 20 else "") + "\n")

    # ── 2. Molecular clock ───────────────────────────────────────────────────
    lines.append(section("2. Molecular clock (root-to-tip divergence)"))
    lines.append(f"| Metric | Value |")
    lines.append(f"|---|---|")
    lines.append(f"| Evolutionary rate | **{slope:.3e}** substitutions/site/year |")
    lines.append(f"| R² | **{r2:.3f}** |")
    lines.append(f"| p-value | {p_val:.2e} |")
    lines.append(f"| Mean root-to-tip divergence | {mean_div:.4f} subs/site |")
    lines.append(f"| Published EV-D68 VP1 rate | ~3–7 × 10⁻³ subs/site/year |")
    lines.append(f"| Rate comparison | {rate_status} |")
    lines.append("")

    r2_interp = interpret_r2(r2)
    lines.append(f"**Clock signal:** {r2_interp}.\n")
    lines.append(
        "![Divergence clock](figures/divergence_clock.png)\n"
        "*Figure 1. Root-to-tip divergence vs. collection date. "
        "The regression slope estimates the evolutionary rate; "
        "R² measures how clock-like the data are.*\n"
    )

    # ── 3. Outliers ──────────────────────────────────────────────────────────
    lines.append(section("3. Clock outliers"))
    if outliers:
        lines.append(
            f"{len(outliers)} sequence(s) deviated more than "
            f"{OUTLIER_SD} SD from the regression line. "
            f"These are highlighted in orange in Figure 1.\n"
        )
        lines.append("| Accession | Year | Date | Divergence | Residual (SD) | Direction |")
        lines.append("|---|---|---|---|---|---|")
        for row in outliers:
            sd_str = row.get("sd_from_line", "—")
            lines.append(
                f"| [{row['accession']}](https://www.ncbi.nlm.nih.gov/nuccore/{row['accession']}) "
                f"| {row['year']} | {row['date']} | {row['divergence']} "
                f"| {sd_str} | {row['direction']} |"
            )
        lines.append("")
        lines.append(
            "> **What to check:** open each accession on NCBI and verify "
            "that the collection date in the record matches the date you would "
            "expect given its position in the tree. Sequences *above* the "
            "regression line evolved faster than average (or are misdated to "
            "an earlier year); sequences *below* evolved slower (or are misdated "
            "to a later year).\n"
        )
    else:
        lines.append(
            f"No outliers detected at the {OUTLIER_SD} SD threshold. "
            "All sequences are consistent with the molecular clock.\n"
        )

    # ── 4. Nucleotide diversity ───────────────────────────────────────────────
    lines.append(section("4. Per-site nucleotide diversity"))
    if H_arr is not None and not np.all(np.isnan(H_arr)):
        mean_H  = float(np.nanmean(H_arr))
        max_H   = float(np.nanmax(H_arr))
        max_pos = int(pos_arr[np.nanargmax(H_arr)])

        # Codon-position means
        pos1_H = H_arr[0::3]
        pos2_H = H_arr[1::3]
        pos3_H = H_arr[2::3]
        m1, m2, m3 = float(np.nanmean(pos1_H)), float(np.nanmean(pos2_H)), float(np.nanmean(pos3_H))

        top5 = top_entropy_sites(H_arr)

        lines.append(f"| Property | Value |")
        lines.append(f"|---|---|")
        lines.append(f"| Alignment positions | {len(H_arr)} |")
        lines.append(f"| Mean entropy | {mean_H:.3f} bits |")
        lines.append(f"| Peak entropy | {max_H:.3f} bits at position {max_pos} |")
        lines.append(f"| Codon pos 1 mean H | {m1:.3f} bits |")
        lines.append(f"| Codon pos 2 mean H | {m2:.3f} bits |")
        lines.append(f"| Codon pos 3 mean H | {m3:.3f} bits (expect highest) |")
        lines.append("")

        lines.append("**Top 5 most variable sites:**\n")
        lines.append("| Rank | Position (1-indexed) | Codon position | Entropy (bits) |")
        lines.append("|---|---|---|---|")
        for rank, (p, h) in enumerate(top5, 1):
            cp = ((p - 1) % 3) + 1
            lines.append(f"| {rank} | {p} | {cp} | {h:.3f} |")
        lines.append("")

        if m3 > m1 and m3 > m2:
            codon_note = (
                "✅ 3rd codon positions are the most variable, as expected "
                "for a protein-coding gene (synonymous mutations are tolerated "
                "more readily than non-synonymous ones)."
            )
        else:
            codon_note = (
                "⚠ 3rd codon positions are *not* the most variable. "
                "This may indicate strong positive selection at specific sites, "
                "an alignment frame offset, or non-coding sequence in the alignment."
            )
        lines.append(f"{codon_note}\n")
        lines.append(
            "![Diversity plot](figures/diversity_plot.png)\n"
            "*Figure 2. Per-site Shannon entropy across the VP1 alignment. "
            "High-entropy sites (top) correspond to rapidly evolving, "
            "potentially immune-selected positions.*\n"
        )
    else:
        lines.append(
            "Diversity data not found. Run `analyze_diversity.py` to generate it.\n"
        )

    # ── 5. Tree ───────────────────────────────────────────────────────────────
    lines.append(section("5. Phylogenetic tree"))
    lines.append(
        "![Phylogenetic tree](figures/phylogenetic_tree.png)\n"
        "*Figure 3. Neighbour-joining tree of EV-D68 VP1 sequences, midpoint-rooted. "
        "Tip colours indicate collection year (blue = earliest, red = most recent).*\n"
    )
    lines.append(
        "The tree was built using the neighbour-joining algorithm "
        "(pairwise identity distance matrix, BioPython `DistanceCalculator`). "
        "Branch lengths are in substitutions per site. "
        "Midpoint rooting was applied as no outgroup sequence was included.\n"
    )
    lines.append(
        "> **Limitation:** Neighbour-joining with an identity distance matrix "
        "is a fast approximation. For a publication, replace this step with "
        "maximum-likelihood inference (e.g. IQ-TREE with a GTR+Γ substitution "
        "model) for more accurate branch lengths and better bootstrapconfidence estimates.\n"
    )

    # ── 6. Methods summary ───────────────────────────────────────────────────
    lines.append(section("6. Methods summary"))
    lines.append(
        "| Step | Tool | Version |"
        "\n|---|---|---|"
        "\n| Sequence retrieval | NCBI Entrez (Biopython `Bio.Entrez`) | Biopython ≥ 1.79 |"
        "\n| Multiple sequence alignment | MAFFT `--auto` | any recent |"
        "\n| Phylogenetic inference | Neighbour-joining (Biopython `DistanceTreeConstructor`) | — |"
        "\n| Root-to-tip analysis | Custom Python (scipy linear regression) | — |"
        "\n| Entropy calculation | Custom Python (Shannon formula) | — |"
        "\n| Visualisation | Matplotlib | ≥ 3.6 |"
        "\n"
    )

    # ── 7. Suggested next steps ───────────────────────────────────────────────
    lines.append(section("7. Suggested next steps"))
    lines.append(
        "1. **IQ-TREE maximum-likelihood tree.** "
        "Replace the NJ tree with `iqtree2 -s data/ev_d68_vp1_aligned.fasta -m GTR+G -bb 1000` "
        "for bootstrap support values and more accurate branch lengths.\n"
        "2. **TempEst root-to-tip validation.** "
        "Open the Newick tree in [TempEst](http://tree.bio.ed.ac.uk/software/tempest/) "
        "to confirm midpoint rooting and compare R² estimates.\n"
        "3. **Protein translation.** Translate VP1 to amino acids and map "
        "the top high-entropy sites from Step 5 onto the known EV-D68 crystal "
        "structure (PDB: 5YHN) to see whether they fall on the antigenic surface.\n"
        "4. **Expand to other genes.** Repeat the pipeline for the 2A protease or "
        "2C helicase to compare evolutionary rates across the genome.\n"
        "5. **Apollo simulation.** Use van Marle et al.'s "
        "[Apollo simulator](https://github.com/van-marle-lab/apollo) to generate "
        "synthetic genomes with a known evolutionary rate, then test whether "
        "this NJ+clock pipeline recovers that rate accurately.\n"
    )

    # ── 8. References ─────────────────────────────────────────────────────────
    lines.append(section("8. References"))
    lines.append(
        "- Tan et al. (2023). Enterovirus D68 evolution and emergence. "
        "*Lancet Microbe*.\n"
        "- Rambaut et al. (2016). Exploring the temporal structure of "
        "heterochronous sequences using TempEst. "
        "*Virus Evolution* 2(1): vew007.\n"
        "- Katoh & Standley (2013). MAFFT multiple sequence alignment. "
        "*Molecular Biology and Evolution* 30(4): 772–780.\n"
        "- Van Marle et al. (2025). Apollo: GPU-accelerated simulation "
        "of viral evolution within a host. *Nature Communications*.\n"
    )

    # ── write file ────────────────────────────────────────────────────────────
    with open(OUT_REPORT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Wrote {OUT_REPORT}")
    print(f"  Sections: dataset, clock, outliers, diversity, tree, methods, "
          f"next steps, references")


if __name__ == "__main__":
    main()
