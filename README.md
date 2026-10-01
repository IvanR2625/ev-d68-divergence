# EV-D68 VP1 Longitudinal Divergence

A Python pipeline that downloads 10+ years of **Enterovirus D68** VP1 capsid
sequences from NCBI, aligns them with MAFFT, builds a neighbour-joining
phylogenetic tree, and plots **root-to-tip divergence against collection date**
— a standard test for molecular-clock behaviour in an RNA virus.

**Why EV-D68?** Enterovirus D68 caused a large outbreak in North America in
2014, and has recurred roughly every two years since. It is an RNA virus with a
mutation rate ~10,000× faster than DNA organisms, making evolutionary signal
visible over just a few years of sampling. VP1 (Viral Protein 1) is the
outermost capsid protein, used to type enteroviruses, and is the standard gene
for evolutionary studies of this group.

---

## What the pipeline produces

| Output | Description |
|--------|-------------|
| `data/ev_d68_vp1_raw.fasta` | Raw VP1 sequences (one per sample, with date + country in header) |
| `data/metadata.csv` | Accession, collection date, country, VP1 length |
| `data/ev_d68_vp1_aligned.fasta` | Multiple sequence alignment (MAFFT) |
| `output/ev_d68_vp1.nwk` | Midpoint-rooted neighbour-joining tree (Newick) |
| `output/ev_d68_vp1_tree.txt` | ASCII tree preview |
| `output/figures/divergence_clock.png` | Three-panel figure (see below) |

### The figure

**Panel 1 — Molecular clock:** scatter plot of root-to-tip divergence
(substitutions per site) vs. collection date, with a linear regression.
The slope is the evolutionary rate in substitutions/site/year; R² measures how
clock-like the data are. RNA viruses typically reach R² > 0.8 with clean
sampling.

**Panel 2 — Sequences per year:** bar chart showing sampling density across the
study period — important context for interpreting the clock plot.

**Panel 3 — Divergence distribution per year:** box plots showing how
within-year divergence changes over time. A rising median confirms directional
evolution.

---

## Repository structure

```
ev-d68-divergence/
├── scripts/
│   ├── fetch_sequences.py    # Step 1 — NCBI Entrez download
│   ├── align_sequences.py    # Step 2 — MAFFT alignment
│   ├── build_tree.py         # Step 3 — NJ tree (BioPython)
│   ├── plot_divergence.py    # Step 4 — divergence figure
│   └── run_pipeline.py       # Run all steps in order
├── data/                     # Created by step 1
├── output/                   # Created by step 3/4
│   └── figures/
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Setup

### 1. Install Python

Download Python 3.10 or later from https://www.python.org/downloads/ — use
the full installer, not the Windows Store version, and tick "Add Python to
PATH" during installation.

### 2. Install Python packages

```bash
pip install -r requirements.txt
```

This installs: `biopython`, `numpy`, `scipy`, `matplotlib`.

### 3. Install MAFFT

MAFFT is a fast multiple sequence aligner — better than pure Python alignment
for hundreds of sequences.

- **Windows:** download the installer from
  https://mafft.cbrc.jp/alignment/software/windows.html
  and follow the instructions. After installing, add the MAFFT folder to your
  system PATH, **or** open `scripts/align_sequences.py` and set the `MAFFT_EXE`
  variable to the full path of `mafft.bat`.
- **macOS:** `brew install mafft`
- **Linux:** `sudo apt install mafft` or `conda install -c bioconda mafft`

---

## Running the pipeline

### One command (recommended)

```bash
python scripts/run_pipeline.py
```

Steps already completed are skipped. To redo everything from scratch:

```bash
python scripts/run_pipeline.py --force
```

### Step by step

```bash
# Step 1 — download sequences from NCBI (~5–15 min depending on connection)
python scripts/fetch_sequences.py

# Step 2 — align with MAFFT (~1–3 min)
python scripts/align_sequences.py

# Step 3 — build neighbour-joining tree (~30 sec)
python scripts/build_tree.py

# Step 4 — plot divergence (~10 sec)
python scripts/plot_divergence.py
```

### Before step 1

Open `scripts/fetch_sequences.py` and set your email address at the top:

```python
ENTREZ_EMAIL = "your.real.email@example.com"
```

NCBI requires an email for Entrez API access. Without it the search will fail
or return an error.

---

## Configuring the download

All adjustable parameters are at the top of each script. The most commonly
changed ones in `fetch_sequences.py`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `YEAR_START` | 2014 | First year to search |
| `YEAR_END` | 2024 | Last year to search |
| `MAX_PER_YEAR` | 15 | Max sequences per year (random sample if more found) |

15 sequences × 11 years = up to 165 sequences. MAFFT handles this in a few
minutes. If you want more resolution, increase `MAX_PER_YEAR` to 30–50; the
alignment and tree steps will take longer.

---

## How each step works

### Step 1 — fetch_sequences.py

Uses NCBI Entrez (via Biopython's `Bio.Entrez`) to search the Nucleotide
database for complete EV-D68 genomes. For each GenBank record it:

1. Finds the `CDS` feature annotated as `gene="VP1"`.
2. Extracts the VP1 nucleotide sequence.
3. Parses the `collection_date` qualifier from the `source` feature.
4. Saves the sequence with the header `>accession|date|country`.

Sequences without a VP1 annotation or without a collection date are skipped.

### Step 2 — align_sequences.py

Runs `mafft --auto --thread -1` on the raw FASTA. The `--auto` flag lets MAFFT
choose the best algorithm based on sequence count and length (FFT-NS-1 for
large inputs, L-INS-i for smaller ones). After alignment, any column consisting
entirely of gap characters (`-`) is removed.

### Step 3 — build_tree.py

Uses Biopython's `DistanceCalculator("identity")` to compute a pairwise
distance matrix (fraction of differing sites), then constructs a
neighbour-joining tree with `DistanceTreeConstructor().nj()`. The tree is
midpoint-rooted: the root is placed at the midpoint of the longest path between
any two leaves, a standard heuristic when no outgroup is available.

The tree is saved in Newick format (`.nwk`), which can be opened in FigTree,
iTOL, or any other tree viewer.

### Step 4 — plot_divergence.py

For each leaf, the root-to-tip distance is computed by summing branch lengths
from the root. These distances are paired with collection dates from
`metadata.csv`. A linear regression (via `scipy.stats.linregress`) estimates
the evolutionary rate. The three-panel figure is saved as a PNG.

---

## Interpreting the results

| What you see | What it means |
|---|---|
| Points trending upward left-to-right | Virus is accumulating mutations over time (expected) |
| High R² (> 0.8) | Strong molecular clock — sampling is consistent and the virus evolves regularly |
| Low R² (< 0.5) | Clock signal is weak — may need more sequences, or the gene is under strong selection |
| Slope ≈ 4–6 × 10⁻³ subs/site/year | Typical for EV-D68 VP1 (published estimates: ~3–7 × 10⁻³) |
| Gap years with no data | Missing from NCBI — add more sequences or adjust `YEAR_START`/`YEAR_END` |

---

## Limitations

- **NJ tree is a rough approximation.** Maximum-likelihood trees (e.g. IQ-TREE)
  give better branch-length estimates and a more reliable clock slope. To
  upgrade step 3, install IQ-TREE and replace the NJ step with:
  `iqtree2 -s data/ev_d68_vp1_aligned.fasta -m GTR+G -bb 1000 --date`.
- **Midpoint rooting can be wrong.** If EV-D68 clades replaced each other
  (which happened in 2014), the midpoint root may land inside a clade rather
  than at the true ancestor. This affects the intercept of the regression but
  not the slope.
- **VP1 only.** Using the full genome (or at minimum the structural proteins
  VP1–VP4) would give more phylogenetic signal.
- **15 sequences/year is sparse.** Increase `MAX_PER_YEAR` for a more robust
  analysis.

---

## Dependencies

| Package | Version | Use |
|---------|---------|-----|
| [Biopython](https://biopython.org/) | ≥ 1.79 | Entrez, SeqIO, alignment, tree building |
| [NumPy](https://numpy.org/) | ≥ 1.23 | Numerical arrays |
| [SciPy](https://scipy.org/) | ≥ 1.9 | Linear regression |
| [Matplotlib](https://matplotlib.org/) | ≥ 3.6 | Figures |
| [MAFFT](https://mafft.cbrc.jp/) | any recent | Multiple sequence alignment (external) |

---

## Further reading

- Tan et al., 2023. "Enterovirus D68 evolution and emergence."
  *Lancet Microbe*. — Overview of EV-D68 molecular evolution since 2014.
- Rambaut et al., 2016. "Exploring the temporal structure of heterochronous
  sequences using TempEst." *Virus Evolution* 2(1): vew007.
  — The tool that does what this pipeline does, more rigorously. Worth
  running on your Newick tree to compare results.
- Van Marle et al., 2025. "Apollo: GPU-accelerated simulation of viral
  evolution." *Nature Communications*.
  — Open-source simulator for viral evolution; the companion paper to
  transmission-inference benchmarking.

---

## License

MIT — use and adapt freely. If you use this pipeline in a science fair project
or publication, a brief acknowledgement is appreciated.
# ev-d68-divergence
# ev-d68-divergence
