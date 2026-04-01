"""
Supplementary Figure 1: DSI threshold sensitivity analysis
============================================================
Paired DSI threshold sweep — dHL-60 compact vs decompact nuclei.

Panel layout (2 × 2):
  [A] AUROC vs threshold
  [B] Cohen's d_z (paired effect size) vs threshold
  [C] Mean ΔDSI (decompact − compact) vs threshold
  [D] Paired DSI at selected τ = 0.3

File naming convention (in compact/ and decompact/ folders):
    <date>_<Dish>_compact_img.tif   / _mask.tif
    <date>_<Dish>_decompact_img.tif / _mask.tif

Author: Minwoo Kang Ph.D., Hawa Thiam Lab, Stanford University
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_auc_score
import tifffile as tiff

# ============================================================
# SETTINGS
# ============================================================
DATA_DIR      = os.path.join(os.path.dirname(__file__), "data")
COMPACT_DIR   = os.path.join(DATA_DIR, "compact")
DECOMPACT_DIR = os.path.join(DATA_DIR, "decompact")

SAVE_SVG  = False
SAVE_DIR  = os.path.join(os.path.dirname(__file__), "outputs")

THRESHOLD_RANGE  = np.linspace(0.05, 0.95, 91)
HIGHLIGHT_THRESH = 0.3

# Figure style
FONT_FAMILY      = 'Arial'
LABEL_SIZE       = 28
TICK_SIZE        = 28
TITLE_SIZE       = 20
SUPTITLE_SIZE    = 24
LEGEND_SIZE      = 15
LINE_WIDTH       = 2.5
AXES_LINEWIDTH   = 2.0
TICK_WIDTH        = 2.0
MARKER_SIZE      = 60
MEAN_MARKER_SIZE = 11
CAPSIZE          = 5
WSPACE           = 0.42
HSPACE           = 0.42

plt.rcParams['font.family']      = FONT_FAMILY
plt.rcParams['axes.linewidth']   = AXES_LINEWIDTH
plt.rcParams['xtick.major.width'] = TICK_WIDTH
plt.rcParams['ytick.major.width'] = TICK_WIDTH
plt.rcParams['xtick.direction']  = 'out'
plt.rcParams['ytick.direction']  = 'out'

# ============================================================
# Helpers
# ============================================================
def minmax_normalize(vals):
    if vals.size == 0:
        return vals
    lo, hi = float(vals.min()), float(vals.max())
    if hi - lo == 0:
        return np.zeros_like(vals)
    return (vals - lo) / (hi - lo)


def get_normalized_pixels(img_path, mask_path):
    img  = tiff.imread(img_path).astype(np.float64)
    mask = tiff.imread(mask_path) > 0
    vals = img[mask]
    if vals.size == 0:
        return np.array([])
    return minmax_normalize(vals)


def compute_dsi(vals_norm, threshold):
    if vals_norm.size == 0:
        return np.nan
    return float(np.mean(vals_norm > threshold))


def cohen_dz(diff):
    if len(diff) < 2:
        return np.nan
    sd = np.std(diff, ddof=1)
    if sd == 0:
        return np.nan
    return float(np.mean(diff) / sd)


def style_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=TICK_SIZE)


def get_prefix(filename):
    base = os.path.basename(filename)
    for keyword in ('_compact_', '_decompact_'):
        if keyword in base:
            return base.split(keyword)[0]
    return '_'.join(base.split('_')[:2])


# ============================================================
# Build matched pairs
# ============================================================
compact_imgs   = sorted(glob.glob(os.path.join(COMPACT_DIR,   "*_compact_img.tif")))
decompact_imgs = sorted(glob.glob(os.path.join(DECOMPACT_DIR, "*_decompact_img.tif")))

compact_dict   = {get_prefix(p): p for p in compact_imgs}
decompact_dict = {get_prefix(p): p for p in decompact_imgs}

common_prefixes = sorted(set(compact_dict) & set(decompact_dict))
print(f"Compact  images found : {len(compact_imgs)}")
print(f"Decompact images found: {len(decompact_imgs)}")
print(f"Matched pairs         : {len(common_prefixes)}")

if len(common_prefixes) == 0:
    raise RuntimeError("No matched pairs found — check filename patterns or folder paths.")

# ============================================================
# Load and normalize all pairs
# ============================================================
px_compact_list   = []
px_decompact_list = []
valid_prefixes    = []

for prefix in common_prefixes:
    c_img_path = compact_dict[prefix]
    d_img_path = decompact_dict[prefix]
    c_mask_path = c_img_path.replace("_img.tif", "_mask.tif")
    d_mask_path = d_img_path.replace("_img.tif", "_mask.tif")

    missing = [p for p in [c_img_path, c_mask_path, d_img_path, d_mask_path]
               if not os.path.exists(p)]
    if missing:
        print(f"  SKIP {prefix} — missing files: {missing}")
        continue

    px_c = get_normalized_pixels(c_img_path, c_mask_path)
    px_d = get_normalized_pixels(d_img_path, d_mask_path)

    if px_c.size == 0 or px_d.size == 0:
        print(f"  SKIP {prefix} — empty ROI after masking")
        continue

    px_compact_list.append(px_c)
    px_decompact_list.append(px_d)
    valid_prefixes.append(prefix)

n_pairs = len(valid_prefixes)
print(f"\nValid pairs used for analysis: {n_pairs}")

# ============================================================
# Threshold sweep
# ============================================================
auroc_vals     = np.full(len(THRESHOLD_RANGE), np.nan)
cohendz_vals   = np.full(len(THRESHOLD_RANGE), np.nan)
wilcox_pvals   = np.full(len(THRESHOLD_RANGE), np.nan)
delta_dsi_mean = np.full(len(THRESHOLD_RANGE), np.nan)

for i, tau in enumerate(THRESHOLD_RANGE):
    dsi_c = np.array([compute_dsi(px, tau) for px in px_compact_list])
    dsi_d = np.array([compute_dsi(px, tau) for px in px_decompact_list])

    delta = dsi_d - dsi_c
    delta_dsi_mean[i] = float(np.mean(delta))

    y_true   = np.concatenate([np.ones(n_pairs), np.zeros(n_pairs)])
    y_scores = np.concatenate([dsi_d, dsi_c])
    try:
        auroc_vals[i] = roc_auc_score(y_true, y_scores)
    except Exception:
        auroc_vals[i] = np.nan

    try:
        _, wilcox_pvals[i] = stats.wilcoxon(dsi_d, dsi_c, alternative='greater')
    except Exception:
        wilcox_pvals[i] = np.nan

    cohendz_vals[i] = cohen_dz(delta)

idx_03 = int(np.argmin(np.abs(THRESHOLD_RANGE - HIGHLIGHT_THRESH)))

# ============================================================
# Summary
# ============================================================
best_auroc_idx    = int(np.nanargmax(auroc_vals))
best_auroc_thresh = THRESHOLD_RANGE[best_auroc_idx]
best_auroc_val    = auroc_vals[best_auroc_idx]

best_dz_idx    = int(np.nanargmax(cohendz_vals))
best_dz_thresh = THRESHOLD_RANGE[best_dz_idx]
best_dz_val    = cohendz_vals[best_dz_idx]

print(f"\n  n pairs               : {n_pairs}")
print(f"  Max AUROC   τ         : {best_auroc_thresh:.2f}  (AUROC={best_auroc_val:.3f})")
print(f"  Max Cohen d_z τ       : {best_dz_thresh:.2f}  (d_z={best_dz_val:.3f})")
print(f"  At selected τ = {HIGHLIGHT_THRESH}:")
print(f"    AUROC     = {auroc_vals[idx_03]:.3f}")
print(f"    Cohen d_z = {cohendz_vals[idx_03]:.3f}")
print(f"    Wilcoxon p= {wilcox_pvals[idx_03]:.2e}")

# ============================================================
# Figure: 2 × 2 panels
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(18, 14))
axes = axes.flatten()

dsi_c_03 = np.array([compute_dsi(px, HIGHLIGHT_THRESH) for px in px_compact_list])
dsi_d_03 = np.array([compute_dsi(px, HIGHLIGHT_THRESH) for px in px_decompact_list])

# --- Panel A: AUROC vs threshold ---
ax = axes[0]
ax.plot(THRESHOLD_RANGE, auroc_vals, color='#333333', linewidth=LINE_WIDTH)
ax.axvline(HIGHLIGHT_THRESH, color='black', linestyle='--', linewidth=2.0,
           label=rf'Selected $\tau$={HIGHLIGHT_THRESH:.1f} (AUC={auroc_vals[idx_03]:.3f})')
ax.axvline(best_auroc_thresh, color='#E55C5C', linestyle=':', linewidth=2.0,
           label=rf'Max AUROC $\tau$={best_auroc_thresh:.2f} (AUC={best_auroc_val:.3f})')
ax.axhline(0.5, color='gray', linewidth=1.0, linestyle=':',
           label='Random classifier baseline (AUROC=0.5)')
ax.set_xlabel(r"DSI threshold ($\tau$)", fontsize=LABEL_SIZE)
ax.set_ylabel("AUROC", fontsize=LABEL_SIZE)
ax.set_title(f"Discriminability\n(n={n_pairs} pairs)", fontsize=TITLE_SIZE)
ax.legend(fontsize=LEGEND_SIZE, frameon=False)
ax.set_xlim(0.05, 0.95); ax.set_ylim(0.4, 1.05)
style_axis(ax)

# --- Panel B: Cohen's d_z vs threshold ---
ax = axes[1]
ax.plot(THRESHOLD_RANGE, cohendz_vals, color='#7B4EA0', linewidth=LINE_WIDTH)
ax.axvline(HIGHLIGHT_THRESH, color='black', linestyle='--', linewidth=2.0,
           label=rf"Selected $\tau$={HIGHLIGHT_THRESH:.1f} (d$_z$={cohendz_vals[idx_03]:.2f})")
ax.axvline(best_dz_thresh, color='#E55C5C', linestyle=':', linewidth=2.0,
           label=rf"Max d$_z$ $\tau$={best_dz_thresh:.2f} (d$_z$={best_dz_val:.2f})")
for d_ref, lbl in [(0.5, 'medium'), (0.8, 'large'), (1.2, 'very large')]:
    ax.axhline(d_ref, color='gray', linewidth=0.8, linestyle=':', alpha=0.7)
    ax.text(0.94, d_ref + 0.03, lbl, transform=ax.get_yaxis_transform(),
            ha='right', fontsize=LEGEND_SIZE, color='gray')
ax.set_xlabel(r"DSI threshold ($\tau$)", fontsize=LABEL_SIZE)
ax.set_ylabel("Cohen's d$_z$ (paired)", fontsize=LABEL_SIZE)
ax.set_title("Paired Effect Size\nvs Threshold", fontsize=TITLE_SIZE)
ax.legend(loc='lower left', fontsize=LEGEND_SIZE, frameon=False)
ax.set_xlim(0.05, 0.95)
style_axis(ax)

# --- Panel C: Mean ΔDSI vs threshold ---
ax = axes[2]
best_delta_idx    = int(np.nanargmax(delta_dsi_mean))
best_delta_thresh = THRESHOLD_RANGE[best_delta_idx]
best_delta_val    = delta_dsi_mean[best_delta_idx]

ax.plot(THRESHOLD_RANGE, delta_dsi_mean, color='#2E8B57', linewidth=LINE_WIDTH)
ax.axvline(HIGHLIGHT_THRESH, color='black', linestyle='--', linewidth=2.0,
           label=rf'Selected $\tau$={HIGHLIGHT_THRESH:.1f} (ΔDSI={delta_dsi_mean[idx_03]:.3f})')
ax.axvline(best_delta_thresh, color='#E55C5C', linestyle=':', linewidth=2.0,
           label=rf'Max ΔDSI $\tau$={best_delta_thresh:.2f} (ΔDSI={best_delta_val:.3f})')
ax.set_xlabel(r"DSI threshold ($\tau$)", fontsize=LABEL_SIZE)
ax.set_ylabel("Mean ΔDSI\n(late − early)", fontsize=LABEL_SIZE)
ax.set_title("Mean Paired Difference\nvs Threshold", fontsize=TITLE_SIZE)
ax.legend(fontsize=LEGEND_SIZE, frameon=False)
ax.set_xlim(0.05, 0.95)
style_axis(ax)

# --- Panel D: Paired DSI at τ = 0.3 ---
ax = axes[3]
for c_val, d_val in zip(dsi_c_03, dsi_d_03):
    ax.plot([1, 2], [c_val, d_val], color='gray', linewidth=0.8, alpha=0.5, zorder=2)

ax.scatter(np.ones(n_pairs), dsi_c_03,
           color='#4A90D9', s=MARKER_SIZE, zorder=4,
           label='Early chromatin reorg.', alpha=0.8)
ax.scatter(np.full(n_pairs, 2), dsi_d_03,
           color='#E5A020', s=MARKER_SIZE, zorder=4,
           label='Late chromatin reorg.', alpha=0.8)

for x_pos, vals, col in [(1, dsi_c_03, '#4A90D9'), (2, dsi_d_03, '#E5A020')]:
    m   = np.mean(vals)
    sem = np.std(vals, ddof=1) / np.sqrt(len(vals))
    ax.errorbar(x_pos, m, yerr=sem, fmt='D', color=col,
                markersize=MEAN_MARKER_SIZE, linewidth=LINE_WIDTH,
                zorder=5, capsize=CAPSIZE)

y_top = max(np.max(dsi_c_03), np.max(dsi_d_03)) + 0.05
ax.plot([1, 2], [y_top, y_top], 'k-', linewidth=1.5)
p_val = wilcox_pvals[idx_03]
stars = ('***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns')
ax.text(1.5, y_top + 0.01, f"{stars}\n(p={p_val:.2e})", ha='center', fontsize=LEGEND_SIZE)

ax.set_xticks([1, 2])
ax.set_xticklabels(['Early chromatin reorg.', 'Late chromatin reorg.'], fontsize=TICK_SIZE)
ax.set_ylabel(rf"DSI ($\tau$ = {HIGHLIGHT_THRESH:.1f})", fontsize=LABEL_SIZE)
ax.set_title(rf"Paired DSI at $\tau$={HIGHLIGHT_THRESH:.1f}" + f"\n(n={n_pairs} nuclei pairs)",
             fontsize=TITLE_SIZE)
ax.legend(fontsize=LEGEND_SIZE, frameon=False)
style_axis(ax)

plt.suptitle(
    "DSI Threshold Sensitivity — dHL-60 Paired Nuclei (early vs late chromatin reorganization)",
    fontsize=SUPTITLE_SIZE, y=0.98, fontweight='bold')

plt.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.90,
                    wspace=WSPACE, hspace=HSPACE)

# ============================================================
# Save
# ============================================================
if __name__ == "__main__":
    os.makedirs(SAVE_DIR, exist_ok=True)

    if SAVE_SVG:
        svg_path = os.path.join(SAVE_DIR, "SuppFig1_threshold_sweep.svg")
        plt.savefig(svg_path, format='svg', bbox_inches='tight')
        print(f"Saved: {svg_path}")
    else:
        plt.show()

    print("Done.")
