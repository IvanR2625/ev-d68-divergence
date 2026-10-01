"""
run_pipeline.py — run all four steps in order.

Usage:
    python scripts/run_pipeline.py

Each step is skipped if its output already exists, so you can re-run
this script safely after fixing a step without redoing earlier work.
Pass --force to redo everything from scratch.

Steps:
  1. fetch_sequences.py  → data/ev_d68_vp1_raw.fasta
  2. align_sequences.py  → data/ev_d68_vp1_aligned.fasta
  3. build_tree.py       → output/ev_d68_vp1.nwk
  4. plot_divergence.py  → output/figures/divergence_clock.png
"""

import os
import sys
import subprocess
import argparse

STEPS = [
    {
        "script": os.path.join("scripts", "fetch_sequences.py"),
        "output": os.path.join("data",   "ev_d68_vp1_raw.fasta"),
        "label":  "Fetch sequences from NCBI",
    },
    {
        "script": os.path.join("scripts", "align_sequences.py"),
        "output": os.path.join("data",   "ev_d68_vp1_aligned.fasta"),
        "label":  "Align with MAFFT",
    },
    {
        "script": os.path.join("scripts", "build_tree.py"),
        "output": os.path.join("output", "ev_d68_vp1.nwk"),
        "label":  "Build neighbour-joining tree",
    },
    {
        "script": os.path.join("scripts", "plot_divergence.py"),
        "output": os.path.join("output", "figures", "divergence_clock.png"),
        "label":  "Plot root-to-tip divergence",
    },
]


def run_step(step, force=False):
    if not force and os.path.exists(step["output"]):
        print(f"  SKIP  {step['label']}  (output exists: {step['output']})")
        return True

    print(f"  RUN   {step['label']} …")
    result = subprocess.run([sys.executable, step["script"]])
    if result.returncode != 0:
        print(f"  FAIL  {step['label']}  (exit code {result.returncode})")
        return False

    print(f"  DONE  {step['label']}")
    return True


def main():
    parser = argparse.ArgumentParser(description="EV-D68 divergence pipeline")
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run all steps even if outputs exist"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("EV-D68 VP1 longitudinal divergence pipeline")
    print("=" * 60)

    for i, step in enumerate(STEPS, 1):
        print(f"\nStep {i}/{len(STEPS)}")
        ok = run_step(step, force=args.force)
        if not ok:
            sys.exit(f"\nPipeline failed at step {i}. Fix the error above and re-run.")

    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print(f"Main output: {STEPS[-1]['output']}")
    print("=" * 60)


if __name__ == "__main__":
    main()
