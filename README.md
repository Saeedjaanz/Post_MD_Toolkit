# Post-MD Analysis Toolkit

**Post-processing and publication-quality analysis of Molecular Dynamics (MD)
simulation trajectories** — works with output from **GROMACS, Desmond, AMBER,
CHARMM, NAMD, OpenMM** and any other engine that MDAnalysis can read.

A single, configurable script reads a topology + trajectory and generates a
complete set of analyses describing the conformational landscape, collective
motions, and protein–ligand contacts of a system:

**PCA · Free Energy Landscape (FEL) · Gibbs Free Energy Landscape · DCCM ·
Contact-frequency bar chart · Contact map (residue × time)**

The toolkit is **system-agnostic** — point it at any protein, protein–ligand, or
protein–protein trajectory and set a few options at the top of the script.

---

## What it produces

| # | Figure | Analysis | Question it answers |
|---|--------|----------|---------------------|
| 1 | `01_PCA.png` | Principal Component Analysis | What are the dominant collective motions? |
| 2 | `02_FEL_3D_surface.png` | Free Energy Landscape (PC1/PC2) | Which conformational states are populated/stable? |
| 3 | `03_Gibbs_FEL_RMSD_Rg.png` | Gibbs Free Energy Landscape (RMSD/Rg) | Is the fold stable & compact through the run? |
| 4 | `04_DCCM.png` | Dynamic Cross-Correlation Matrix | Which residues move together / against each other? |
| 5 | `05_contact_bar.png` | Contact-frequency bar chart | Which residues bind the ligand most often? |
| 6 | `06_contact_heatmap.png` | Contact map (residue × time) | When during the run is each contact formed? |

A companion table, `contact_freq.csv`, lists the per-residue ligand-contact
frequency. The images in [`example_figures/`](example_figures/) are sample
outputs from a demonstration run.

> **A plain-language, step-by-step explanation of every analysis (PCA, FEL,
> Gibbs FEL, DCCM, contacts, contact bar plot) is in
> [`docs/ANALYSIS_GUIDE.md`](docs/ANALYSIS_GUIDE.md).**

---

## Supported MD engines

The script loads trajectories through **MDAnalysis**, so it accepts the native
output of essentially every major MD package. Typical topology + trajectory pairs:

| Engine | Topology | Trajectory |
|--------|----------|------------|
| **GROMACS** | `.gro`, `.tpr` | `.xtc`, `.trr` |
| **AMBER** | `.prmtop`, `.parm7` | `.nc`, `.mdcrd`, `.dcd` |
| **Desmond** | `.cms` | `.dtr` / `_trj` (or convert to `.dcd`/`.xtc`) |
| **CHARMM / NAMD** | `.psf` | `.dcd` |
| **OpenMM** | `.pdb`, `.psf` | `.dcd`, `.xtc` |

Just set `TOPOLOGY` and `TRAJECTORY` in the **CONFIG** block to your files — no
other change is needed.

> **Tip (Desmond):** if the raw `_trj` directory is awkward to load, convert it to
> DCD/XTC first (e.g. with `trj_convert` / VMD), then point the script at that.

---

## Repository layout

```
Post-MD-Analysis-Toolkit/
├── README.md                  # this file
├── requirements.txt           # Python dependencies
├── .gitignore                 # excludes large MD trajectories/binaries
├── src/
│   └── md_analysis.py         # the analysis script (config block at top)
├── docs/
│   └── ANALYSIS_GUIDE.md      # step-by-step guide to each analysis
└── example_figures/           # sample outputs from a demonstration run
    ├── 01_PCA.png
    ├── 02_FEL_3D_surface.png
    ├── 03_Gibbs_FEL_RMSD_Rg.png
    ├── 04_DCCM.png
    ├── 05_contact_bar.png
    ├── 06_contact_heatmap.png
    └── contact_freq.csv
```

---

## Usage

```bash
# 1. install dependencies
pip install -r requirements.txt

# 2. edit the CONFIG block at the top of src/md_analysis.py:
#      TOPOLOGY        -> your topology file   (.gro/.prmtop/.psf/.cms ...)
#      TRAJECTORY      -> your trajectory file  (.xtc/.nc/.dcd ...)
#      LIGAND_RESNAME  -> residue name of your ligand (e.g. UNL, LIG, MOL)

# 3. run
python src/md_analysis.py
```

All figures are written to a `Figures/` directory at 600 DPI.

### Key configuration knobs (top of `src/md_analysis.py`)

| Setting | Default | Meaning |
|---------|---------|---------|
| `TOPOLOGY` | `step5_production.gro` | Topology / coordinates file |
| `TRAJECTORY` | `step5_production_center.xtc` | Trajectory file |
| `LIGAND_RESNAME` | `UNL` | Residue name of the ligand in the topology |
| `T_KELVIN` | `310.0` | Temperature used for Boltzmann inversion |
| `CONTACT_CUTOFF_A` | `4.0` | Heavy-atom distance (Å) defining a contact |
| `CONTACT_MIN_FREQ` | `0.05` | Minimum contact frequency to plot a residue |
| `FEL_BINS` / `GIBBS_BINS` | `100` | 2-D histogram resolution for the energy landscapes |
| `DPI` | `600` | Output resolution |

> **Protein-only run (no ligand)?** Comment out the contact analyses (Figures 5–6)
> or set `LIGAND_RESNAME` to a residue you don't need — PCA, FEL, Gibbs FEL, and
> DCCM run without a ligand.

---

## Method summary

- **Alignment** — the trajectory is least-squares fitted onto frame 0 using Cα atoms.
- **Single pass** — Cα positions, radius of gyration, and ligand contacts are
  collected in one traversal of the trajectory for efficiency.
- **PCA** — performed on mean-centred Cα coordinates (scikit-learn).
- **Free energies** — obtained by Boltzmann inversion of 2-D histograms,
  ΔG = −k_B·T·ln(P).
- **DCCM** — normalised covariance of Cα displacements (vectorised).
- **Contacts** — a residue is "in contact" in a frame when any heavy atom lies
  within the cutoff of any ligand heavy atom.

Full details for each step are in [`docs/ANALYSIS_GUIDE.md`](docs/ANALYSIS_GUIDE.md).

---

## Citation

If you use this toolkit, please cite the underlying libraries:

- **MDAnalysis** — Michaud-Agrawal *et al.*, *J. Comput. Chem.* (2011); Gowers *et al.* (2016)
- **scikit-learn** — Pedregosa *et al.*, *JMLR* (2011)
- Plus the MD engine used to generate the trajectory (GROMACS, AMBER, Desmond, …).

---

## License

Released under the [MIT License](LICENSE) — free to use, modify, and distribute.
