"""
run_pipeline.py — run all pipeline steps in order.

Usage:
    python scripts/run_pipeline.py             # skip completed steps
    python scripts/run_pipeline.py --force     # redo everything
    python scripts/run_pipeline.py --from 3    # start at step N

Steps:
  1. fetch_sequences.py    → data/ev_d68_vp1_raw.fasta
  2. align_sequences.py    → data/ev_d68_vp1_aligned.fasta
  3. build_tree.py         → output/ev_d68_vp1.nwk
  4. plot_divergence.py    → output/figures/divergence_clock.png
                             output/outliers.csv
  5. analyze_diversity.py  → output/diversity_entropy.csv
                             output/figures/diversity_plot.png
  6. plot_tree.py          → output/figures/phylogenetic_tree.png
  7. generate_report.py    → output/REPORT.md
"""

import os
import sys
import subprocess
import argparse

STEPS = [
    {
        "n":      1,
        "script": os.path.join("scripts", "fetch_sequences.py"),
        "output": os.path.join("data",    "ev_d68_vp1_raw.fasta"),
        "label":  "Fetch sequences from NCBI",
    },
    {
        "n":      2,
        "script": os.path.join("scripts", "align_sequences.py"),
        "output": os.path.join("data",    "ev_d68_vp1_aligned.fasta"),
        "label":  "Align with MAFFT",
    },
    {
        "n":      3,
        "script": os.path.join("scripts", "build_tree.py"),
        "output": os.path.join("output",  "ev_d68_vp1.nwk"),
        "label":  "Build neighbour-joining tree",
    },
    {
        "n":      4,
        "script": os.path.join("scripts", "plot_divergence.py"),
        "output": os.path.join("output",  "figures", "divergence_clock.png"),
        "label":  "Plot root-to-tip divergence + outlier detection",
    },
    {
        "n":      5,
        "script": os.path.join("scripts", "analyze_diversity.py"),
        "output": os.path.join("output",  "figures", "diversity_plot.png"),
        "label":  "Per-site entropy and codon-position diversity",
    },
    {
        "n":      6,
        "script": os.path.join("scripts", "plot_tree.py"),
        "output": os.path.join("output",  "figures", "phylogenetic_tree.png"),
        "label":  "Draw year-coloured phylogenetic tree",
    },
    {
        "n":      7,
        "script": os.path.join("scripts", "generate_report.py"),
        "output": os.path.join("output",  "REPORT.md"),
        "label":  "Generate Markdown results report",
    },
]


def run_step(step, force=False):
    label  = f"Step {step['n']}: {step['label']}"
    output = step["output"]

    if not force and os.path.exists(output):
        print(f"  SKIP  {label}")
        print(f"        (output exists: {output})")
        return True

    print(f"\n  ── {label} ──")
    result = subprocess.run([sys.executable, step["script"]])
    if result.returncode != 0:
        print(f"\n  ✗ FAILED  {label}  (exit code {result.returncode})")
        return False

    print(f"  ✓ DONE   {label}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="EV-D68 VP1 longitudinal divergence pipeline"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run all steps even if outputs already exist",
    )
    parser.add_argument(
        "--from", dest="from_step", type=int, default=1, metavar="N",
        help="Start from step N (1–7)",
    )
    parser.add_argument(
        "--only", type=int, default=None, metavar="N",
        help="Run only step N",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  EV-D68 VP1 longitudinal divergence pipeline")
    print(f"  {len(STEPS)} steps total")
    print("=" * 60)

    for step in STEPS:
        n = step["n"]
        if args.only is not None and n != args.only:
            continue
        if n < args.from_step:
            print(f"  SKIP  Step {n}: {step['label']}  (before --from {args.from_step})")
            continue

        ok = run_step(step, force=args.force)
        if not ok:
            print(
                f"\nPipeline stopped at step {n}.\n"
                "Fix the error above, then re-run with:\n"
                f"  python scripts/run_pipeline.py --from {n}"
            )
            sys.exit(1)

    print("\n" + "=" * 60)
    print("  Pipeline complete.")
    print(f"  Results report: {STEPS[-1]['output']}")
    print(f"  Figures:        output/figures/")
    print("=" * 60)


if __name__ == "__main__":
    main()
