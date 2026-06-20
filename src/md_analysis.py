#!/usr/bin/env python3
# =============================================================================
# MD Analyses: PCA, FEL, Gibbs FEL, DCCM, Contact Heatmap
# PCA figures restyled with QuickPCA cosmetics (TheVisualHub/QuickPCA)
# System : generic protein–ligand MD (set per run in the CONFIG block below)
# Author : Saeed Ishaq
# Date   : 2026-05-04
# =============================================================================

# =============================================================================
# ██████╗ ██████╗ ███╗   ██╗███████╗██╗ ██████╗
# ██╔════╝██╔═══██╗████╗  ██║██╔════╝██║██╔════╝
# ██║     ██║   ██║██╔██╗ ██║█████╗  ██║██║  ███╗
# ██║     ██║   ██║██║╚██╗██║██╔══╝  ██║██║   ██║
# ╚██████╗╚██████╔╝██║ ╚████║██║     ██║╚██████╔╝
#  ╚═════╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝     ╚═╝ ╚═════╝
# =============================================================================

# ─── INPUT FILES ──────────────────────────────────────────────────────────────
TOPOLOGY       = "step5_production.gro"
TRAJECTORY     = "step5_production_center.xtc"
OUTPUT_DIR     = "Figures"

# ─── SYSTEM DEFINITION ────────────────────────────────────────────────────────
# Set LIGAND_RESNAME to your ligand's residue name for THIS system (e.g. UNL,
# LIG, MOL). It is validated against the topology at load time. Chain/Cα
# composition is auto-detected at runtime — nothing here is system-specific.
LIGAND_RESNAME     = "UNL"
CA_SELECT          = "protein and name CA"
PROTEIN_HEAVY_SEL  = "protein and not name H*"
LIGAND_HEAVY_SEL   = f"resname {LIGAND_RESNAME} and not name H*"

# ─── THERMODYNAMIC PARAMETERS ─────────────────────────────────────────────────
T_KELVIN  = 310.0
KB_KJ_MOL = 8.314e-3
KBT       = KB_KJ_MOL * T_KELVIN    # 2.4942 kJ/mol

# ─── PCA SETTINGS ─────────────────────────────────────────────────────────────
N_PCA_COMPONENTS  = 10
PCA_ALIGN_FRAME   = 0
PCA_SCATTER_ALPHA = 0.55

# ─── FEL / GIBBS FEL SETTINGS ─────────────────────────────────────────────────
FEL_BINS       = 100
FEL_MAX_ENERGY = 20.0
GIBBS_BINS     = 100
FEL_SMOOTH_SIG = 1.2    # gaussian_filter sigma for FEL smoothing

# ─── CONTACT HEATMAP SETTINGS ─────────────────────────────────────────────────
CONTACT_CUTOFF_A = 4.0
CONTACT_MIN_FREQ = 0.05

# ─── FIGURE / PUBLICATION QUALITY ─────────────────────────────────────────────
DPI             = 600
FIG_FORMAT      = ["png"]
FONT_FAMILY     = "Liberation Serif"
FONT_SIZE       = 13
TICK_SIZE       = 11
LEGEND_SIZE     = 10
TITLE_SIZE      = 14
CBAR_LABEL_SIZE = 12
LINE_WIDTH      = 1.8
MARKER_SIZE     = 16
SPINE_WIDTH     = 1.2

# ─── COLORMAPS ────────────────────────────────────────────────────────────────
CMAP_FEL      = "viridis"
CMAP_DCCM     = "RdBu_r"
CMAP_CONTACT  = "YlOrRd"
CMAP_PCA_TIME = "viridis"

# ─── COLOURS — all black, no grey text ────────────────────────────────────────
C_BLACK   = "#000000"
C_LINE1   = "#1F5C99"
C_LINE2   = "#B22222"
C_THRESH  = "#8B0000"

# ─── QuickPCA cosmetic palette (PCA figures) ──────────────────────────────────
CMAP_PCA_FEL_Q = "RdYlBu_r"     # QuickPCA free-energy landscape colormap
Q_START        = "lime"         # trajectory start marker
Q_END          = "red"          # trajectory end marker
Q_PC1          = "teal"         # PC1 distribution
Q_PC2          = "darkorange"   # PC2 distribution
Q_BAR          = "steelblue"    # scree bars
Q_BAR_EDGE     = "navy"         # scree bar edges
Q_CUM          = "coral"        # cumulative-variance line

# =============================================================================
# END CONFIG
# =============================================================================

import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.ndimage import gaussian_filter
from MDAnalysis.analysis import align, rms, distances

warnings.filterwarnings("ignore")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─── GLOBAL MATPLOTLIB STYLE ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"         : FONT_FAMILY,
    "font.size"           : FONT_SIZE,
    "text.color"          : C_BLACK,
    "axes.labelcolor"     : C_BLACK,
    "axes.edgecolor"      : C_BLACK,
    "axes.titlecolor"     : C_BLACK,
    "xtick.color"         : C_BLACK,
    "ytick.color"         : C_BLACK,
    "xtick.labelcolor"    : C_BLACK,
    "ytick.labelcolor"    : C_BLACK,
    "axes.titlesize"      : TITLE_SIZE,
    "axes.labelsize"      : FONT_SIZE,
    "xtick.labelsize"     : TICK_SIZE,
    "ytick.labelsize"     : TICK_SIZE,
    "legend.fontsize"     : LEGEND_SIZE,
    "legend.edgecolor"    : C_BLACK,
    "legend.framealpha"   : 0.93,
    "axes.linewidth"      : SPINE_WIDTH,
    "xtick.major.width"   : SPINE_WIDTH,
    "ytick.major.width"   : SPINE_WIDTH,
    "xtick.minor.width"   : 0.9,
    "ytick.minor.width"   : 0.9,
    "xtick.major.size"    : 5,
    "ytick.major.size"    : 5,
    "xtick.minor.size"    : 3,
    "ytick.minor.size"    : 3,
    "xtick.direction"     : "in",
    "ytick.direction"     : "in",
    "axes.spines.top"     : True,
    "axes.spines.right"   : True,
    "axes.facecolor"      : "#FAFAFA",
    "figure.facecolor"    : "white",
    "figure.dpi"          : DPI,
    "savefig.dpi"         : DPI,
    "savefig.bbox"        : "tight",
    "savefig.facecolor"   : "white",
    "savefig.transparent" : False,
    "axes.titleweight"    : "bold",
    "axes.labelweight"    : "bold",
})


def _style_ax(ax, minor=True):
    for sp in ax.spines.values():
        sp.set_linewidth(SPINE_WIDTH)
        sp.set_edgecolor(C_BLACK)
    ax.tick_params(which="both", direction="in",
                   top=True, right=True,
                   colors=C_BLACK, labelcolor=C_BLACK)
    if minor:
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))


def _cbar(fig, im, ax, label, orient="vertical", frac=0.046):
    cb = fig.colorbar(im, ax=ax, orientation=orient, pad=0.025, fraction=frac)
    cb.set_label(label, fontsize=CBAR_LABEL_SIZE, color=C_BLACK, fontweight="bold")
    cb.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 1)
    cb.outline.set_edgecolor(C_BLACK)
    cb.outline.set_linewidth(SPINE_WIDTH)
    return cb


def _fel_grid(G, bins):
    """histogram → normalize → NaN-mask zeros → Boltzmann inversion → smooth."""
    h, xe, ye = np.histogram2d(*bins, bins=FEL_BINS if G == "fel" else GIBBS_BINS)
    hn = h / h.max()
    hn[hn == 0] = np.nan
    G_raw = -KBT * np.log(hn)
    G_raw[G_raw > FEL_MAX_ENERGY] = np.nan
    G_sm  = gaussian_filter(np.nan_to_num(G_raw, nan=FEL_MAX_ENERGY), sigma=FEL_SMOOTH_SIG)
    G_sm[np.isnan(G_raw)] = np.nan
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    mi, mj = np.unravel_index(np.nanargmin(G_raw), G_raw.shape)
    return G_sm, xc, yc, mi, mj


def _plot_fel(fig, ax, G_sm, xc, yc, mi, mj, xlabel, ylabel, title):
    im = ax.pcolormesh(xc, yc, G_sm.T, cmap=CMAP_FEL,
                       shading="nearest", rasterized=True, zorder=1)
    im.set_clim(0, FEL_MAX_ENERGY)
    _cbar(fig, im, ax, "ΔG (kJ mol⁻¹)")
    try:
        ax.contour(xc, yc, G_sm.T, levels=[2, 4, 6, 8, 10, 14],
                   colors=C_BLACK, linewidths=0.75, alpha=0.5, zorder=2)
    except Exception:
        pass
    ax.plot(xc[mi], yc[mj], "*", ms=MARKER_SIZE,
            color="#FFD700", mew=1.4, mec=C_BLACK, zorder=5,
            label="Global minimum  (ΔG = 0)")
    ax.legend(loc="upper right", frameon=True,
              edgecolor=C_BLACK, fontsize=LEGEND_SIZE)
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.set_title(title, fontweight="bold", pad=10)
    _style_ax(ax)


def save_fig(fig, name):
    for fmt in FIG_FORMAT:
        path = os.path.join(OUTPUT_DIR, f"{name}.{fmt}")
        fig.savefig(path, dpi=DPI)
        print(f"  Saved → {path}")
    plt.close(fig)


# =============================================================================
# LOAD UNIVERSE & SANITY CHECKS
# =============================================================================
print("=" * 70)
print("Loading universe …")
u   = mda.Universe(TOPOLOGY, TRAJECTORY)
ref = mda.Universe(TOPOLOGY, TRAJECTORY)

assert LIGAND_RESNAME in set(u.residues.resnames), \
    f"Ligand '{LIGAND_RESNAME}' not found."

ca_atoms   = u.select_atoms(CA_SELECT)
lig_heavy  = u.select_atoms(LIGAND_HEAVY_SEL)
prot_heavy = u.select_atoms(PROTEIN_HEAVY_SEL)

assert len(ca_atoms)   > 0, "Empty Cα selection."
assert len(lig_heavy)  > 0, "Empty ligand heavy-atom selection."
assert len(prot_heavy) > 0, "Empty protein heavy-atom selection."

n_ca      = ca_atoms.n_atoms
n_frames  = u.trajectory.n_frames
ca_resids = ca_atoms.resids
time_ns   = np.array([u.trajectory[i].time / 1000.0 for i in range(n_frames)])

# Auto-detect chain composition from the topology (segids). Some formats
# (e.g. .gro) do not store chain IDs, in which case everything is one segment.
_ca_segids = ca_atoms.segids
_seen      = list(dict.fromkeys(_ca_segids))          # unique, order-preserving
chain_ca   = [(sid, int((_ca_segids == sid).sum())) for sid in _seen]
if len(chain_ca) > 1:
    chain_str = ", ".join(f"{sid or '?'}:{n}" for sid, n in chain_ca)
else:
    chain_str = "single segment (chain IDs not stored in topology)"

print(f"  Frames   : {n_frames}")
print(f"  Time span: 0 – {time_ns[-1]:.0f} ns")
print(f"  Cα       : {n_ca}  ({chain_str})")
print(f"  Ligand   : {LIGAND_RESNAME}  ({len(lig_heavy)} heavy atoms)")

try:
    masses = u.select_atoms("protein").masses
    print(f"  Masses   : {'BAD (zeros/NaN)' if np.any(masses==0)|np.any(np.isnan(masses)) else 'OK'}")
except:
    print("  Masses   : WARNING – cannot read")


# =============================================================================
# STEP 1 — ALIGN TRAJECTORY (Cα, frame 0)
# =============================================================================
print("\n[1/5] Aligning trajectory …")
ref.trajectory[PCA_ALIGN_FRAME]
align.AlignTraj(u, ref, select=CA_SELECT, in_memory=True).run()
print("  Done.")


# =============================================================================
# STEP 2 — SINGLE-PASS DATA COLLECTION
# Collect Cα positions, Rg, and contact counts in ONE trajectory traversal.
# =============================================================================
print("\n[2/5] Single-pass data collection …")

protein_residues = prot_heavy.residues
n_prot_res       = len(protein_residues)
res_heavy_list   = [res.atoms.select_atoms("not name H*") for res in protein_residues]
prot_sel         = u.select_atoms("protein")

pos_nm          = np.zeros((n_frames, n_ca * 3), dtype=np.float32)   # PCA (nm)
pos_ca          = np.zeros((n_frames, n_ca, 3),  dtype=np.float32)   # DCCM (Å)
rg_nm           = np.zeros(n_frames, dtype=np.float64)
contact_cnt     = np.zeros(n_prot_res, dtype=np.int32)
contact_matrix  = np.zeros((n_frames, n_prot_res), dtype=np.bool_)   # per-frame binary

for i, ts in enumerate(u.trajectory):
    ca_pos = ca_atoms.positions
    pos_nm[i]  = ca_pos.flatten() / 10.0
    pos_ca[i]  = ca_pos

    try:
        rg_nm[i] = prot_sel.radius_of_gyration() / 10.0
    except Exception:
        pos_p = prot_sel.positions
        rg_nm[i] = np.sqrt(np.mean(np.sum((pos_p - pos_p.mean(0))**2, axis=1))) / 10.0

    lig_pos = lig_heavy.positions
    for j, ratoms in enumerate(res_heavy_list):
        if distances.distance_array(lig_pos, ratoms.positions).min() <= CONTACT_CUTOFF_A:
            contact_cnt[j] += 1
            contact_matrix[i, j] = True

    if (i + 1) % 500 == 0 or i == n_frames - 1:
        print(f"    frame {i+1}/{n_frames}")

contact_freq = contact_cnt / n_frames
print("  Collection complete.")


# =============================================================================
# STEP 3 — RMSD  (reuse aligned trajectory; MDAnalysis RMSD uses cached coords)
# =============================================================================
print("\n[3/5] Computing RMSD …")
ref.trajectory[PCA_ALIGN_FRAME]
rmsd_run = rms.RMSD(ca_atoms, ref.select_atoms(CA_SELECT),
                    select=CA_SELECT).run()
rmsd_nm  = rmsd_run.results.rmsd[:, 2] / 10.0    # Å → nm
print(f"  RMSD: {rmsd_nm.min():.3f}–{rmsd_nm.max():.3f} nm  "
      f"Rg: {rg_nm.min():.3f}–{rg_nm.max():.3f} nm")


# =============================================================================
# STEP 4 — PCA
# =============================================================================
print("\n[4/5] PCA …")
pca_model  = PCA(n_components=N_PCA_COMPONENTS)
pca_model.fit(pos_nm)
pc_proj    = pca_model.transform(pos_nm)
evr        = pca_model.explained_variance_ratio_ * 100.0
pc1, pc2   = pc_proj[:, 0], pc_proj[:, 1]
print(f"  PC1: {evr[0]:.1f}%   PC2: {evr[1]:.1f}%   Top-{N_PCA_COMPONENTS}: {evr.cumsum()[-1]:.1f}%")


# =============================================================================
# STEP 5 — DCCM  (vectorised — no Python nested loop)
# =============================================================================
print("\n[5/5] DCCM …")
mean_ca  = pos_ca.mean(axis=0)                               # (N_CA, 3)
delta    = (pos_ca - mean_ca).reshape(n_frames, n_ca * 3).astype(np.float64)
cov_3N   = (delta.T @ delta) / n_frames                     # (3N, 3N)

# Vectorised block-trace extraction: reshape cov to (N_CA, 3, N_CA, 3)
# then sum diagonal of each 3×3 block → (N_CA, N_CA)
cov_4d   = cov_3N.reshape(n_ca, 3, n_ca, 3)
cov_atom = np.einsum("iajb,ab->ij", cov_4d,
                     np.eye(3, dtype=np.float64))            # (N_CA, N_CA)

var_d    = np.sqrt(np.diag(cov_atom))
denom    = np.outer(var_d, var_d)
denom[denom == 0] = np.nan
dccm     = cov_atom / denom
assert np.allclose(dccm, dccm.T, atol=1e-5, equal_nan=True), "DCCM not symmetric!"
print(f"  Range: {np.nanmin(dccm):.3f} – {np.nanmax(dccm):.3f}")


# =============================================================================
# FREE MEMORY before plotting
# =============================================================================
del pos_nm, pos_ca, delta, cov_3N, cov_4d, cov_atom
import gc; gc.collect()


# =============================================================================
# FIGURES
# =============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — PCA  (Panel A style: solid viridis FEL surface, kcal/mol)
# ──────────────────────────────────────────────────────────────────────────────
KJ_TO_KCAL = 0.239006

G_fel, xc_f, yc_f, mi_f, mj_f = _fel_grid("fel", (pc1, pc2))

# --- kept (kcal/mol, solid-filled) for the 3D surface in Figure 3 ---
G_kcal       = G_fel * KJ_TO_KCAL
MAX_KCAL     = FEL_MAX_ENERGY * KJ_TO_KCAL
G_kcal_solid = np.where(np.isnan(G_kcal), MAX_KCAL, G_kcal)
G_kcal_solid = gaussian_filter(G_kcal_solid, sigma=1.5)

# --- 2D FEL: density histogram → Boltzmann inversion → shift min=0 ---
# (rebuilt to match QuickPCA's colorbar: levels to 97th percentile + extend="max")
H_q, xe_q, ye_q = np.histogram2d(pc1, pc2, bins=FEL_BINS, density=True)
H_q = gaussian_filter(H_q, sigma=FEL_SMOOTH_SIG)
with np.errstate(divide="ignore", invalid="ignore"):
    F_q = np.where(H_q > 0, -KBT * np.log(H_q), np.nan)
F_q -= np.nanmin(F_q)                       # global minimum → 0

xc_q = 0.5 * (xe_q[:-1] + xe_q[1:])
yc_q = 0.5 * (ye_q[:-1] + ye_q[1:])
XX_q, YY_q = np.meshgrid(xc_q, yc_q, indexing="ij")

F_max_q  = np.nanpercentile(F_q, 97)        # top of colour scale (~ like ref's 42)
levels_q = np.linspace(0, F_max_q, 30)
F_plot   = np.where(np.isnan(F_q), F_max_q, F_q)   # solid high-energy background

fig, ax = plt.subplots(figsize=(8.0, 6.8))

cf = ax.contourf(XX_q, YY_q, F_plot, levels=levels_q,
                 cmap=CMAP_PCA_FEL_Q, extend="max")
try:
    ax.contour(XX_q, YY_q, F_plot, levels=levels_q[::4],
               colors="k", linewidths=0.4, alpha=0.5, zorder=2)
except Exception:
    pass

# faint white trajectory line through the projected frames + start/end stars
ax.plot(pc1, pc2, color="white", lw=0.25, alpha=0.3, zorder=3)
ax.scatter(pc1[0],  pc2[0],  c=Q_START, s=130, marker="*",
           edgecolors="k", linewidths=0.7, zorder=5, label="Start")
ax.scatter(pc1[-1], pc2[-1], c=Q_END,   s=130, marker="*",
           edgecolors="k", linewidths=0.7, zorder=5, label="End")

cb = fig.colorbar(cf, ax=ax, pad=0.025, fraction=0.046)
cb.set_label("Free Energy (kJ mol⁻¹)", fontsize=CBAR_LABEL_SIZE,
             color=C_BLACK, fontweight="bold")
_step = max(1, int(round(F_max_q / 6)))     # ~6 evenly-spaced integer ticks
cb.set_ticks(range(0, int(F_max_q) + _step, _step))
cb.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 1)
cb.outline.set_edgecolor(C_BLACK); cb.outline.set_linewidth(SPINE_WIDTH)

ax.set_xlabel(f"PC1 ({evr[0]:.1f}%)", fontweight="bold")
ax.set_ylabel(f"PC2 ({evr[1]:.1f}%)", fontweight="bold")
ax.set_title("Principal Component Analysis", fontweight="bold", pad=10)
# Clamp axes to the filled FEL extent — removes the ~5% margin gap at the borders
ax.set_xlim(xc_q.min(), xc_q.max())
ax.set_ylim(yc_q.min(), yc_q.max())
ax.grid(color="white", alpha=0.15, linestyle="--", linewidth=0.5)
ax.legend(loc="upper right", frameon=True, edgecolor=C_BLACK, fontsize=LEGEND_SIZE)
_style_ax(ax, minor=False)
fig.tight_layout()
save_fig(fig, "01_PCA")


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — FEL  single panel 3D solid surface  (viridis, kcal/mol)
# ──────────────────────────────────────────────────────────────────────────────
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

XX_f, YY_f = np.meshgrid(xc_f, yc_f, indexing="ij")

fig = plt.figure(figsize=(9, 7.5), facecolor="white")
ax3f = fig.add_subplot(111, projection="3d")

surf_f = ax3f.plot_surface(XX_f, YY_f, G_kcal_solid,
                            cmap=CMAP_PCA_FEL_Q, vmin=0, vmax=MAX_KCAL,
                            linewidth=0, antialiased=True,
                            rstride=1, cstride=1, rasterized=True)

cb_3f = fig.colorbar(surf_f, ax=ax3f, pad=0.12, fraction=0.030, shrink=0.7)
cb_3f.set_label("Free Energy (kcal mol⁻¹)", fontsize=CBAR_LABEL_SIZE,
                color=C_BLACK, fontweight="bold")
cb_3f.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 2)
cb_3f.outline.set_edgecolor(C_BLACK); cb_3f.outline.set_linewidth(SPINE_WIDTH)

ax3f.set_xlabel("PC1 (nm)", fontsize=TICK_SIZE, labelpad=8)
ax3f.set_ylabel("PC2 (nm)", fontsize=TICK_SIZE, labelpad=8)
ax3f.set_zlabel("ΔG (kcal mol⁻¹)",   fontsize=TICK_SIZE, labelpad=8)
ax3f.set_title("Free Energy Landscape (FEL)", fontweight="bold", pad=14,
               fontsize=TITLE_SIZE)
ax3f.tick_params(labelsize=TICK_SIZE - 2, colors=C_BLACK)
ax3f.xaxis.pane.set_facecolor("white"); ax3f.xaxis.pane.set_edgecolor(C_BLACK)
ax3f.yaxis.pane.set_facecolor("white"); ax3f.yaxis.pane.set_edgecolor(C_BLACK)
ax3f.zaxis.pane.set_facecolor("white"); ax3f.zaxis.pane.set_edgecolor(C_BLACK)
ax3f.xaxis.pane.fill = True
ax3f.yaxis.pane.fill = True
ax3f.zaxis.pane.fill = True
ax3f.view_init(elev=32, azim=-55)

fig.tight_layout()
save_fig(fig, "02_FEL_3D_surface")


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — Gibbs FEL  two-panel  (2D pcolormesh + 3D surface)
# Colormap: jet  (blue=low, red=high);  unsampled bins filled with max → red bg
# ──────────────────────────────────────────────────────────────────────────────
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

G_gib, xc_g, yc_g, mi_g, mj_g = _fel_grid("gib", (rmsd_nm, rg_nm))

# Fill NaN → max energy (solid red bg in jet); then smooth for continuous surface
G_gib_fill   = np.where(np.isnan(G_gib), FEL_MAX_ENERGY, G_gib)
G_gib_smooth = gaussian_filter(G_gib_fill, sigma=1.2)
XX_g, YY_g   = np.meshgrid(xc_g, yc_g, indexing="ij")

fig = plt.figure(figsize=(15, 6.4), facecolor="white")

# — left: 2D pcolormesh (gouraud for smooth look) —
ax2d = fig.add_subplot(121)
_style_ax(ax2d)
im_g2 = ax2d.pcolormesh(xc_g, yc_g, G_gib_smooth.T, cmap="jet",
                         shading="gouraud", rasterized=True)
im_g2.set_clim(0, FEL_MAX_ENERGY)
cb_g = fig.colorbar(im_g2, ax=ax2d, pad=0.025, fraction=0.046)
cb_g.set_label("Gibbs Free Energy (kJ mol⁻¹)", fontsize=CBAR_LABEL_SIZE,
               color=C_BLACK, fontweight="bold")
cb_g.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 1)
cb_g.outline.set_edgecolor(C_BLACK); cb_g.outline.set_linewidth(SPINE_WIDTH)
ax2d.set_xlabel("RMSD (nm)", fontweight="bold")
ax2d.set_ylabel("R$_g$ (nm)", fontweight="bold")
ax2d.set_title("Gibbs Free Energy Landscape", fontweight="bold", pad=10)

# — right: 3D solid surface —
ax3d = fig.add_subplot(122, projection="3d")
surf = ax3d.plot_surface(XX_g, YY_g, G_gib_smooth,
                          cmap="jet", vmin=0, vmax=FEL_MAX_ENERGY,
                          linewidth=0, antialiased=True,
                          rstride=1, cstride=1, rasterized=True)
cb3d = fig.colorbar(surf, ax=ax3d, pad=0.12, fraction=0.030, shrink=0.7)
cb3d.set_label("Gibbs Free Energy (kJ mol⁻¹)", fontsize=CBAR_LABEL_SIZE - 1,
               color=C_BLACK, fontweight="bold")
cb3d.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 2)
cb3d.outline.set_edgecolor(C_BLACK); cb3d.outline.set_linewidth(SPINE_WIDTH)
ax3d.set_xlabel("RMSD (nm)",      fontsize=TICK_SIZE, labelpad=7)
ax3d.set_ylabel("R$_g$ (nm)",     fontsize=TICK_SIZE, labelpad=7)
ax3d.set_zlabel("ΔG (kJ mol⁻¹)", fontsize=TICK_SIZE, labelpad=7)
ax3d.set_title("Gibbs Free Energy Landscape", fontweight="bold", pad=12)
ax3d.tick_params(labelsize=TICK_SIZE - 2, colors=C_BLACK)
# Solid white panes so the base is visible
ax3d.xaxis.pane.set_facecolor("white"); ax3d.xaxis.pane.set_edgecolor(C_BLACK)
ax3d.yaxis.pane.set_facecolor("white"); ax3d.yaxis.pane.set_edgecolor(C_BLACK)
ax3d.zaxis.pane.set_facecolor("white"); ax3d.zaxis.pane.set_edgecolor(C_BLACK)
ax3d.xaxis.pane.fill = True
ax3d.yaxis.pane.fill = True
ax3d.zaxis.pane.fill = True
ax3d.view_init(elev=30, azim=-60)

fig.tight_layout()
save_fig(fig, "03_Gibbs_FEL_RMSD_Rg")


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — DCCM
# ──────────────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9.5, 8.5))
_style_ax(ax, minor=False)

im_d = ax.imshow(dccm, cmap=CMAP_DCCM, vmin=-1, vmax=1,
                 aspect="auto", interpolation="bilinear",
                 origin="lower", zorder=1)
cb_d = _cbar(fig, im_d, ax, "Cross-correlation coefficient  ($C_{ij}$)")
cb_d.set_ticks([-1, -0.5, 0, 0.5, 1])

step = max(1, n_ca // 10)
tpos = list(range(0, n_ca, step))
tlab = [str(ca_resids[i]) for i in tpos]
ax.set_xticks(tpos); ax.set_xticklabels(tlab, rotation=45, ha="right",
                                          fontsize=TICK_SIZE, color=C_BLACK)
ax.set_yticks(tpos); ax.set_yticklabels(tlab, fontsize=TICK_SIZE, color=C_BLACK)
ax.set_xlabel("Residue number", fontweight="bold")
ax.set_ylabel("Residue number", fontweight="bold")
ax.set_title("Dynamic Cross-Correlation Matrix (DCCM)", fontweight="bold", pad=14)
fig.tight_layout()
save_fig(fig, "04_DCCM")


# ──────────────────────────────────────────────────────────────────────────────
# CONTACT DATA
# ──────────────────────────────────────────────────────────────────────────────
res_resids = np.array([r.resid   for r in protein_residues])
res_names  = np.array([r.resname for r in protein_residues])

df_all = pd.DataFrame({"resid": res_resids, "resname": res_names,
                        "freq": contact_freq})
df_all.to_csv(os.path.join(OUTPUT_DIR, "contact_freq.csv"), index=False)
print(f"\nContacts — non-zero: {(contact_freq > 0).sum()}"
      f"   significant (≥30%): {(contact_freq >= 0.30).sum()}")

mask    = contact_freq >= CONTACT_MIN_FREQ
df_plot = df_all[mask].copy().reset_index(drop=True)
cmap_c  = plt.cm.get_cmap(CMAP_CONTACT)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 5 — Contact Bar Chart
# ──────────────────────────────────────────────────────────────────────────────
if len(df_plot) > 0:
    x_pos    = np.arange(len(df_plot))
    x_labels = [f"{r.resname}{r.resid}" for _, r in df_plot.iterrows()]
    bar_cols = cmap_c(df_plot["freq"].values)

    fig_w = max(9, len(df_plot) * 0.38)
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))
    _style_ax(ax)

    ax.bar(x_pos, df_plot["freq"].values, color=bar_cols,
           edgecolor=C_BLACK, linewidth=0.6, width=0.72, zorder=2)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.55, alpha=0.4,
                  color=C_BLACK, zorder=0)
    ax.set_axisbelow(True)

    ax.axhline(0.30, color=C_THRESH, lw=1.6, ls="--", zorder=3,
               label="30% significance threshold (ProLIF)")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=55, ha="right",
                       fontsize=9, color=C_BLACK)
    ax.set_ylim(0, 1.14)
    ax.set_ylabel(f"Contact frequency  (cutoff = {CONTACT_CUTOFF_A} Å)",
                  fontweight="bold")
    ax.set_xlabel("Protein residue", fontweight="bold")
    ax.set_title(
        f"Protein–Ligand Contact Frequency  (LIG  |  ≥{CONTACT_MIN_FREQ*100:.0f}% frames shown)",
        fontweight="bold", pad=10)

    sm = ScalarMappable(cmap=CMAP_CONTACT, norm=Normalize(0, 1))
    sm.set_array([])
    _cbar(fig, sm, ax, "Contact frequency", frac=0.025)
    ax.legend(loc="upper left", frameon=True,
              edgecolor=C_BLACK, fontsize=LEGEND_SIZE)
    fig.tight_layout()
    save_fig(fig, "05_contact_bar")


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 6 — Contact Heatmap  (residues × frames)
# Y-axis : three-letter residue code + residue number  (e.g. GLU226)
# X-axis : simulation frames
# Colour : red–yellow  (light yellow = no contact, deep red = contact)
# ──────────────────────────────────────────────────────────────────────────────
from matplotlib.colors import LinearSegmentedColormap

# Custom red-yellow colormap: white → light yellow → orange → deep red
CMAP_RY = LinearSegmentedColormap.from_list(
    "ry", ["#FFFFFF", "#FFF59D", "#FF8F00", "#C62828"], N=256
)

if len(df_plot) > 0:
    # Indices of significant residues in the full residue array
    sig_idx   = np.where(mask)[0]                         # positions in protein_residues
    ct_sub    = contact_matrix[:, sig_idx].T.astype(float)  # (n_sig_res × n_frames)
    n_sig     = len(sig_idx)

    # Y-axis labels: three-letter code + resid  (e.g. "GLU 226")
    y_labels  = [f"{r.resname} {r.resid}" for r in [protein_residues[k] for k in sig_idx]]

    # Figure height scales with number of residues (0.32 in each)
    fig_h = max(4.0, n_sig * 0.32 + 1.5)
    fig_w = max(12.0, n_frames / 400.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im_c = ax.imshow(ct_sub, cmap=CMAP_RY, vmin=0, vmax=1,
                     aspect="auto", interpolation="nearest",
                     origin="upper")

    # Colorbar
    cb = fig.colorbar(im_c, ax=ax, orientation="vertical",
                      pad=0.015, fraction=0.018)
    cb.set_label("Contact  (1 = present)", fontsize=CBAR_LABEL_SIZE,
                 color=C_BLACK, fontweight="bold")
    cb.set_ticks([0, 1])
    cb.set_ticklabels(["Absent", "Present"])
    cb.ax.tick_params(labelcolor=C_BLACK, color=C_BLACK, labelsize=TICK_SIZE - 1)
    cb.outline.set_edgecolor(C_BLACK); cb.outline.set_linewidth(SPINE_WIDTH)

    # Y-axis — residue labels
    ax.set_yticks(np.arange(n_sig))
    ax.set_yticklabels(y_labels, fontsize=8.5, color=C_BLACK,
                       fontfamily=FONT_FAMILY)

    # X-axis — frame numbers with time (ns) labels
    n_xticks  = 10
    xt_idx    = np.linspace(0, n_frames - 1, n_xticks, dtype=int)
    xt_labels = [f"{time_ns[k]:.0f}" for k in xt_idx]
    ax.set_xticks(xt_idx)
    ax.set_xticklabels(xt_labels, fontsize=TICK_SIZE, color=C_BLACK)

    ax.set_xlabel("Simulation time (ns)", fontweight="bold", color=C_BLACK)
    ax.set_ylabel("Residue", fontweight="bold", color=C_BLACK)
    ax.set_title(
        "Protein–Ligand Contact Map",
        fontweight="bold", pad=10, color=C_BLACK)

    for sp in ax.spines.values():
        sp.set_edgecolor(C_BLACK); sp.set_linewidth(SPINE_WIDTH)
    ax.tick_params(which="both", direction="in",
                   colors=C_BLACK, labelcolor=C_BLACK,
                   top=True, right=True)

    fig.tight_layout()
    save_fig(fig, "06_contact_heatmap")


# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("ALL ANALYSES COMPLETE")
print(f"Output: {os.path.abspath(OUTPUT_DIR)}/")
print("-" * 70)
for f in sorted(os.listdir(OUTPUT_DIR)):
    print(f"  {f}")
print("=" * 70)
