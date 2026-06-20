# Analysis Guide — what each step does

A short, to-the-point explanation of every analysis in `src/md_analysis.py`,
in the order the script runs them. Each section says **what** it is, **how** the
script computes it, and **how to read** the resulting figure.

---

## 0. Setup (alignment & single-pass data collection)

Before any analysis, the trajectory is **aligned** so that overall translation and
rotation of the protein are removed — only *internal* motion remains.

- The protein is least-squares fitted onto **frame 0** using its **Cα atoms**
  (`align.AlignTraj`, `in_memory=True`).
- Then the script walks the trajectory **once** and, for every frame, stores:
  - Cα coordinates (for PCA and DCCM),
  - radius of gyration `Rg` (compactness),
  - which protein residues are within `4 Å` of the ligand (contacts).

Doing everything in one pass keeps a large trajectory fast to process.

> **Why align?** PCA and DCCM measure how atoms move *relative to the rest of the
> protein*. If you don't remove whole-molecule tumbling, those motions swamp the
> real internal dynamics.

---

## 1. PCA — Principal Component Analysis → `01_PCA.png`

**What:** PCA finds the handful of collective motions ("principal components",
PCs) that explain most of the protein's movement. Instead of 3·N noisy
coordinates, you describe the dynamics with 2–3 numbers per frame.

**How the script does it:**
- Input = mean-centred Cα coordinates (in nm) for every frame.
- `sklearn.decomposition.PCA` is fit; the top 10 components are kept.
- Each frame is projected onto **PC1** and **PC2**.
- `explained_variance_ratio_` reports how much motion each PC captures.

**How to read the figure:**
- Axes are **PC1** and **PC2**; the % on each axis is the variance it explains.
- The coloured surface is a **free-energy map** of the projection
  (blue/low = frequently visited & stable; red/high = rarely visited).
- The **green star** = first frame, **red star** = last frame, so you can see
  how far the conformation drifted over the run.
- One deep basin = a single stable state; multiple basins = the protein
  hops between conformations.

---

## 2. FEL — Free Energy Landscape (PC1 vs PC2) → `02_FEL_3D_surface.png`

**What:** the same PC1/PC2 space as Figure 1, but drawn as a **3-D energy
surface** so the depth of each conformational basin is obvious.

**How the script does it (Boltzmann inversion):**
1. Build a 2-D histogram of (PC1, PC2) → probability `P` of each region.
2. Convert probability to free energy with
   **ΔG = −k_B·T·ln(P/P_max)**, with `T = 310 K`.
3. Unsampled / very-high-energy bins are capped at `FEL_MAX_ENERGY` and the
   surface is smoothed (Gaussian filter) for a clean plot.
4. The global minimum is shifted to ΔG = 0.

**How to read the figure:**
- **Valleys (low ΔG, blue)** = stable, populated conformations.
- **Peaks (high ΔG, red)** = unfavourable, rarely sampled states.
- A single broad/deep valley → the complex settled into one stable structure.
- Units here are **kcal/mol** (the 2-D version in Fig 1 uses kJ/mol).

---

## 3. Gibbs FEL — RMSD vs Rg → `03_Gibbs_FEL_RMSD_Rg.png`

**What:** a free-energy landscape using two intuitive physical coordinates instead
of abstract PCs:
- **RMSD** — how far the structure has moved from the starting structure
  (conformational drift).
- **Rg (radius of gyration)** — how compact/folded the protein is.

**How the script does it:** identical Boltzmann inversion as Fig 2, but the 2-D
histogram is over (RMSD, Rg). Shown as a **2-D map (left)** and a **3-D surface
(right)** with the `jet` colormap.

**How to read the figure:**
- The energy minimum tells you the **preferred (RMSD, Rg)** of the complex.
- **Low, tight RMSD + stable Rg** = the protein stays folded and near its start
  → a stable, well-behaved simulation.
- Several minima = the protein samples distinct folded/expanded states.

---

## 4. DCCM — Dynamic Cross-Correlation Matrix → `04_DCCM.png`

**What:** a residue-by-residue map of **correlated motion**. It tells you which
parts of the protein move *together* and which move in *opposite* directions.

**How the script does it:**
- Compute the covariance of Cα displacements about their mean position.
- Normalise to a correlation coefficient
  **C_ij = ⟨Δr_i·Δr_j⟩ / (σ_i·σ_j)**, ranging from −1 to +1
  (done vectorised with `einsum`, then checked for symmetry).

**How to read the figure (RdBu_r colormap):**
- **+1 (deep blue)** = residues i and j move **in the same direction**
  (correlated).
- **−1 (deep red)** = they move in **opposite directions** (anti-correlated).
- **0 (white)** = independent motion.
- The diagonal is always +1 (every residue with itself).
- Off-diagonal blue blocks reveal **domains** moving as rigid units; red blocks
  reveal **hinge / opposing** motions (often functionally important).

---

## 5. Contact bar chart → `05_contact_bar.png`

**What:** ranks the protein residues by how often they touch the ligand across the
whole trajectory — i.e. the **binding-site fingerprint**.

**How the script does it:**
- During the single pass, a residue counts as "in contact" in a frame if **any of
  its heavy atoms** is within `CONTACT_CUTOFF_A = 4 Å` of any ligand heavy atom.
- **Contact frequency** = (frames in contact) / (total frames).
- Only residues with frequency ≥ `CONTACT_MIN_FREQ` (5 %) are plotted.

**How to read the figure:**
- Each bar = one residue (labelled `RESNAME+number`, e.g. `GLU226`).
- Bar height & colour (YlOrRd) = contact frequency (taller/darker = more
  persistent contact).
- The dashed line at **0.30** marks a **30 % significance threshold**: residues
  above it are considered key, stable interaction partners of the ligand.

---

## 6. Contact heatmap / map → `06_contact_heatmap.png`

**What:** the *time-resolved* version of Figure 5 — it shows **when** during the
simulation each contact is present, not just how often.

**How the script does it:** the per-frame binary contact matrix (residue × frame)
is plotted for the significant residues only.

**How to read the figure:**
- **Y-axis** = significant residues (`RESNAME + number`).
- **X-axis** = simulation time (ns).
- **Colour** — light = no contact, deep red = in contact.
- A solid horizontal red band = a residue that stays bound the **entire** run
  (a stable anchor). Broken/intermittent bands = transient contacts that form and
  break.

---

## Quick glossary

| Term | Meaning |
|------|---------|
| **Cα** | The alpha-carbon of each residue; used as a backbone tracer. |
| **RMSD** | Root-mean-square deviation — average atomic displacement from a reference. |
| **Rg** | Radius of gyration — a measure of how compact the protein is. |
| **PC** | Principal component — a dominant collective motion from PCA. |
| **ΔG / FEL** | Free energy / free-energy landscape; low = stable, high = unstable. |
| **Boltzmann inversion** | Turning a probability distribution into free energy via ΔG = −k_B·T·ln(P). |
| **DCCM** | Dynamic cross-correlation matrix of residue motions (−1 to +1). |
| **Contact** | Heavy-atom of a residue within 4 Å of a ligand heavy atom. |
| **k_B·T** | Thermal energy; here 2.4942 kJ/mol at 310 K. |
