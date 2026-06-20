# SERPINE1–Cajanone — MD Trajectory Analysis

Post-processing and publication-quality analysis of a molecular-dynamics (MD)
simulation of the **SERPINE1** protein (a liver-cancer / PAI-1 target) in complex
with the natural-product ligand **Cajanone**.

A single script reads a GROMACS production trajectory and generates a complete set
of figures describing the conformational landscape, collective motions, and
protein–ligand contacts of the complex.

---

## What the script produces

| # | Figure | Analysis | Question it answers |
|---|--------|----------|---------------------|
| 1 | `01_PCA.png` | Principal Component Analysis | What are the dominant collective motions? |
| 2 | `02_FEL_3D_surface.png` | Free Energy Landscape (PC1/PC2) | Which conformational states are populated/stable? |
| 3 | `03_Gibbs_FEL_RMSD_Rg.png` | Gibbs Free Energy Landscape (RMSD/Rg) | Is the fold stable & compact through the run? |
| 4 | `04_DCCM.png` | Dynamic Cross-Correlation Matrix | Which residues move together / against each other? |
| 5 | `05_contact_bar.png` | Contact-frequency bar chart | Which residues bind the ligand most often? |
| 6 | `06_contact_heatmap.png` | Contact map (residue × time) | When during the run is each contact formed? |

A companion table, `figures/contact_freq.csv`, lists the per-residue
ligand-contact frequency.

> **A plain-language, step-by-step explanation of every analysis (PCA, FEL,
> Gibbs FEL, DCCM, contacts, contact bar plot) is in
> [`docs/ANALYSIS_GUIDE.md`](docs/ANALYSIS_GUIDE.md).**

---

## Repository layout

```
SERPINE1-Cajanone-MD-Analysis/
├── README.md                  # this file
├── requirements.txt           # Python dependencies
├── .gitignore                 # excludes large MD trajectories/binaries
├── src/
│   └── md_analysis.py         # the analysis script (config block at top)
├── docs/
│   └── ANALYSIS_GUIDE.md      # step-by-step guide to each analysis
└── figures/
    ├── 01_PCA.png
    ├── 02_FEL_3D_surface.png
    ├── 03_Gibbs_FEL_RMSD_Rg.png
    ├── 04_DCCM.png
    ├── 05_contact_bar.png
    ├── 06_contact_heatmap.png
    └── contact_freq.csv
```

---

## Inputs (not committed)

The script needs two files from the GROMACS run. They are **not** included in this
repo because MD trajectories are large (GB-scale); see `.gitignore`.

| Variable | File | Description |
|----------|------|-------------|
| `TOPOLOGY` | `step5_production.gro` | Topology / coordinates of the production system |
| `TRAJECTORY` | `step5_production_center.xtc` | Centred production trajectory |

Place both files next to the script (or edit the paths in the **CONFIG** block) and
run.

---

## Usage

```bash
# 1. install dependencies
pip install -r requirements.txt

# 2. set the system in the CONFIG block of src/md_analysis.py
#    most importantly LIGAND_RESNAME (default "UNL")

# 3. run from the directory holding the .gro and .xtc
python src/md_analysis.py
```

All figures are written to a `Figures/` directory at 600 DPI.

### Key configuration knobs (top of `src/md_analysis.py`)

| Setting | Default | Meaning |
|---------|---------|---------|
| `LIGAND_RESNAME` | `UNL` | Residue name of the ligand in the topology |
| `T_KELVIN` | `310.0` | Temperature used for Boltzmann inversion (body temp) |
| `CONTACT_CUTOFF_A` | `4.0` | Heavy-atom distance (Å) defining a contact |
| `CONTACT_MIN_FREQ` | `0.05` | Minimum contact frequency to plot a residue |
| `FEL_BINS` / `GIBBS_BINS` | `100` | 2-D histogram resolution for the energy landscapes |
| `DPI` | `600` | Output resolution |

---

## Method summary

- **Alignment** — the trajectory is least-squares fitted onto frame 0 using Cα atoms.
- **Single pass** — Cα positions, radius of gyration, and ligand contacts are
  collected in one traversal of the trajectory for efficiency.
- **PCA** — performed on mean-centred Cα coordinates (scikit-learn).
- **Free energies** — obtained by Boltzmann inversion of 2-D histograms,
  ΔG = −k_B·T·ln(P), with T = 310 K.
- **DCCM** — normalised covariance of Cα displacements (vectorised).
- **Contacts** — a residue is "in contact" in a frame when any heavy atom lies
  within 4 Å of any ligand heavy atom.

Full details for each step are in [`docs/ANALYSIS_GUIDE.md`](docs/ANALYSIS_GUIDE.md).

---

## Citation

If you use this analysis, please cite the underlying tools:

- **MDAnalysis** — Michaud-Agrawal *et al.*, *J. Comput. Chem.* (2011); Gowers *et al.* (2016)
- **scikit-learn** — Pedregosa *et al.*, *JMLR* (2011)
- **GROMACS** — Abraham *et al.*, *SoftwareX* (2015)

Author of the analysis script: **Saeed Ishaq**.
