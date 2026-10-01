# EV-D68 VP1 Longitudinal Divergence

A Python pipeline that downloads 10+ years of **Enterovirus D68** VP1 capsid
sequences from NCBI, aligns them with MAFFT, builds a phylogenetic tree, and
produces a suite of evolutionary analyses — all from publicly available data,
no lab access required.

**Why EV-D68?** Enterovirus D68 caused a large outbreak in North America in
2014, and has recurred roughly every two years since. As an RNA virus its
mutation rate is ~10,000× faster than DNA organisms, making evolutionary signal
visible over just a few years of sampling. VP1 (Viral Protein 1) is the
outermost capsid protein — it contacts antibodies, undergoes immune selection,
and is the standard gene for evolutionary studies of enteroviruses.

**Authors:** Ivan Raizada  
**Affiliation:** Webber Academy, Calgary, AB

---

## Outputs

| File | What it shows |
|------|---------------|
| `output/figures/divergence_clock.png` | Root-to-tip divergence vs. date — tests the molecular clock; slope = evolutionary rate |
| `output/figures/diversity_plot.png` | Per-site Shannon entropy across VP1 — reveals immune-selected vs. conserved regions |
| `output/figures/phylogenetic_tree.png` | NJ tree with tips coloured by collection year |
| `output/outliers.csv` | Sequences that deviate > 2 SD from the clock regression |
| `output/diversity_entropy.csv` | Entropy at every aligned position |
| `output/REPORT.md` | Auto-generated Markdown report with all key statistics |

---

## Pipeline overview

```
NCBI Nucleotide
      │  fetch_sequences.py
      ▼
data/ev_d68_vp1_raw.fasta
      │  align_sequences.py  (MAFFT)
      ▼
data/ev_d68_vp1_aligned.fasta ──────────┐
      │  build_tree.py                  │ analyze_diversity.py
      ▼                                 ▼
output/ev_d68_vp1.nwk        output/diversity_entropy.csv
      │                                 │
      ├─ plot_divergence.py             │ diversity_plot.png
      │    divergence_clock.png         │
      │    outliers.csv                 │
      │                                 │
      └─ plot_tree.py                   │
           phylogenetic_tree.png        │
                    │                   │
                    └───────────────────┘
                               │  generate_report.py
                               ▼
                        output/REPORT.md
```

---

## Scripts

| Script | Step | Input | Output |
|--------|------|-------|--------|
| `fetch_sequences.py` | 1 | NCBI Entrez API | `data/ev_d68_vp1_raw.fasta`, `data/metadata.csv` |
| `align_sequences.py` | 2 | raw FASTA | `data/ev_d68_vp1_aligned.fasta` |
| `build_tree.py` | 3 | aligned FASTA | `output/ev_d68_vp1.nwk`, `output/ev_d68_vp1_tree.txt` |
| `plot_divergence.py` | 4 | tree + metadata | `divergence_clock.png`, `outliers.csv` |
| `analyze_diversity.py` | 5 | aligned FASTA | `diversity_entropy.csv`, `diversity_plot.png` |
| `plot_tree.py` | 6 | tree + metadata | `phylogenetic_tree.png` |
| `generate_report.py` | 7 | all outputs | `REPORT.md` |
| `utils.py` | — | shared helpers | imported by steps 4–7 |
| `run_pipeline.py` | — | orchestrator | runs steps 1–7 in order |

---

## Setup

### 1. Install Python

Download Python **3.10 or later** from https://www.python.org/downloads/ — use the
full installer (not the Windows Store version) and tick **"Add Python to PATH"**
during installation.

### 2. Install Python packages

```bash
pip install -r requirements.txt
```

Installs: `biopython`, `numpy`, `scipy`, `matplotlib`.

### 3. Install MAFFT

MAFFT is the external aligner used in step 2.

| Platform | Command |
|----------|---------|
| **Windows** | Download the installer from https://mafft.cbrc.jp/alignment/software/windows.html, install, then add the install folder to PATH — **or** set `MAFFT_EXE` in `align_sequences.py` to the full path of `mafft.bat`. |
| **macOS** | `brew install mafft` |
| **Linux** | `sudo apt install mafft` or `conda install -c bioconda mafft` |

### 4. Set your NCBI email

Open `scripts/fetch_sequences.py` and replace the placeholder:

```python
ENTREZ_EMAIL = "your.real.email@example.com"
```

NCBI requires an email for Entrez API access. The script exits immediately if
you leave the placeholder in place.

---

## Running the pipeline

### One command

```bash
python scripts/run_pipeline.py
```

Steps already completed are skipped. Re-run options:

```bash
python scripts/run_pipeline.py --force        # redo everything
python scripts/run_pipeline.py --from 4       # restart from step 4
python scripts/run_pipeline.py --only 5       # run just step 5
```

### Step by step

```bash
python scripts/fetch_sequences.py     # ~5–15 min (network)
python scripts/align_sequences.py     # ~1–3 min
python scripts/build_tree.py          # ~30 sec
python scripts/plot_divergence.py     # ~10 sec
python scripts/analyze_diversity.py   # ~10 sec
python scripts/plot_tree.py           # ~20 sec
python scripts/generate_report.py     # ~5 sec
```

---

## How each step works

### Step 1 — fetch_sequences.py

Uses NCBI Entrez (`Bio.Entrez`) to search for complete EV-D68 genomes
year by year. For each GenBank record it:

1. Finds the `CDS` feature annotated with `gene="VP1"`.
2. Extracts the VP1 nucleotide sequence.
3. Parses the `collection_date` qualifier from the `source` feature.
4. Saves the sequence as `>accession|YYYY-MM-DD|Country`.

Sequences without a VP1 annotation or a parseable collection date are dropped.
Up to `MAX_PER_YEAR` sequences are kept per year (default 15; random sample
with a fixed seed for reproducibility).

### Step 2 — align_sequences.py

Runs `mafft --auto --thread -1` on the raw FASTA. `--auto` lets MAFFT choose
the best algorithm for the input size. After alignment, any column that is all
gaps is removed.

### Step 3 — build_tree.py

Computes a pairwise identity distance matrix
(`Bio.Phylo.TreeConstruction.DistanceCalculator`) and builds a
neighbour-joining tree. The tree is midpoint-rooted using
`tree.root_at_midpoint()`. Output is saved as Newick plus an ASCII preview.

### Step 4 — plot_divergence.py

For each tip, the root-to-tip distance is computed by summing branch lengths
from the root. Dates are parsed from sequence headers. A linear regression
(`scipy.stats.linregress`) gives the evolutionary rate (slope) and clock
quality (R²). Sequences whose residual exceeds `OUTLIER_SD × std(residuals)`
are flagged as clock outliers and highlighted in orange. A 95% confidence band
is drawn around the regression line.

**Outlier interpretation:** sequences above the regression line evolved faster
than average (or are misdated to an earlier year); below means slower (or
misdated to a later year). Outliers are saved to `output/outliers.csv` with
NCBI links for manual verification.

### Step 5 — analyze_diversity.py

For each column in the aligned FASTA, nucleotide frequencies are computed
(gaps and ambiguous bases ignored). Shannon entropy is:

```
H = -Σ p_i × log₂(p_i)    for nucleotides with p_i > 0
```

H = 0 (all sequences identical) to H = 2 (uniform A/T/G/C). A 30-bp sliding
window average highlights regional patterns. Entropy is also split by codon
position (1, 2, 3) — 3rd positions should have the highest entropy because
synonymous changes are tolerated without affecting protein function.

### Step 6 — plot_tree.py

Reads the Newick tree and draws it with `Bio.Phylo.draw()`. Tip labels are
shortened to `YYYY_accession` and coloured with the `coolwarm` colormap scaled
to the year range. A colorbar legend maps colour to year.

### Step 7 — generate_report.py

Reads all available output files and writes `output/REPORT.md` — a complete
Markdown report with a dataset table, clock statistics, outlier table
(with NCBI links), diversity summary, tree description, methods, and suggested
next steps.

---

## Configuring the download

All adjustable parameters are at the top of each script. The most useful ones:

**`fetch_sequences.py`**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ENTREZ_EMAIL` | placeholder | **Required.** Your NCBI email address. |
| `YEAR_START` | 2014 | First year to search |
| `YEAR_END` | 2024 | Last year to search |
| `MAX_PER_YEAR` | 15 | Max sequences per year |

**`plot_divergence.py`**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `OUTLIER_SD` | 2.0 | SD threshold for clock outlier detection |

**`analyze_diversity.py`**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WINDOW_SIZE` | 30 | Sliding window width (nt) for smoothed entropy |

---

## Interpreting the results

### Molecular clock (Step 4)

| What you see | What it means |
|---|---|
| Points trending upward | Virus is accumulating mutations over time ✓ |
| R² ≥ 0.85 | Strong clock — consistent evolutionary rate |
| R² 0.70–0.84 | Reasonable — more sequences per year would help |
| R² < 0.50 | Weak — check for sampling gaps or misdated records |
| Slope ≈ 3–7 × 10⁻³ | In the published range for EV-D68 VP1 |
| Orange diamonds | Clock outliers — verify their dates on NCBI |

### Per-site diversity (Step 5)

| What you see | What it means |
|---|---|
| High-entropy peaks | Rapidly evolving, likely immune-selected residues |
| Low-entropy baseline | Conserved sites — mutations here disrupt protein function |
| 3rd codon positions most variable | Expected for protein-coding genes |
| 3rd positions NOT most variable | Unusual — strong positive selection, or check alignment frame |

### Tree (Step 6)

| What you see | What it means |
|---|---|
| Same-colour tips spread throughout tree | Good temporal mixing — sequences from each year throughout |
| Long terminal branches on recent (red) tips | Fast-evolving lineages or recent clade emergence |
| 'Broom' pattern at tips | Population expansion |
| Two distinct clades | Multiple co-circulating EV-D68 lineages (known for this virus) |

---

## Limitations and next steps

- **NJ tree is an approximation.** For publication, use IQ-TREE with a GTR+Γ
  model:
  ```bash
  iqtree2 -s data/ev_d68_vp1_aligned.fasta -m GTR+G -bb 1000
  ```
- **Midpoint rooting can be wrong** when multiple clades co-circulate, as they
  do for EV-D68. An outgroup (e.g. EV-D70 VP1) would give a more reliable root.
- **15 sequences/year is sparse.** Increase `MAX_PER_YEAR` for a more robust
  rate estimate.
- **VP1 only.** Including all four structural proteins (VP1–VP4) would give
  more phylogenetic signal.
- **SMOTE / class balance** does not apply here — all sequences are included.

**Suggested further analyses** (also listed in `output/REPORT.md`):

1. Validate the clock in [TempEst](http://tree.bio.ed.ac.uk/software/tempest/)
   using the generated Newick file.
2. Map the top high-entropy sites from Step 5 onto the EV-D68 crystal structure
   (PDB: 5YHN) to confirm they fall on antigenic surface loops.
3. Run [Apollo](https://github.com/van-marle-lab/apollo) (van Marle et al., 2025)
   to generate synthetic genomes with a known rate and test whether this pipeline
   recovers it — a direct benchmarking exercise.
4. Extend to other EV-D68 genes (2A protease, 2C helicase) to compare rates
   across the genome.

---

## Dependencies

| Package | Version | Use |
|---------|---------|-----|
| [Biopython](https://biopython.org/) | ≥ 1.79 | Entrez, SeqIO, alignment, tree |
| [NumPy](https://numpy.org/) | ≥ 1.23 | Numerical arrays |
| [SciPy](https://scipy.org/) | ≥ 1.9 | Linear regression, confidence intervals |
| [Matplotlib](https://matplotlib.org/) | ≥ 3.6 | All figures |
| [MAFFT](https://mafft.cbrc.jp/) | any recent | Multiple sequence alignment (external binary) |

---

## Key references

- Tan et al. (2023). Enterovirus D68 evolution and emergence. *Lancet Microbe*.
- Rambaut et al. (2016). Exploring the temporal structure of heterochronous
  sequences using TempEst. *Virus Evolution* 2(1): vew007. — the tool this
  pipeline's clock analysis is based on.
- Katoh & Standley (2013). MAFFT multiple sequence alignment software.
  *Molecular Biology and Evolution* 30(4): 772–780.
- Van Marle et al. (2025). Apollo: GPU-accelerated simulation of viral evolution
  within a host. *Nature Communications*. — open-source viral evolution simulator;
  a natural companion to this pipeline.

---

## License

MIT — use and adapt freely. If you use this pipeline in a science fair project
or publication, a brief acknowledgement is appreciated.
