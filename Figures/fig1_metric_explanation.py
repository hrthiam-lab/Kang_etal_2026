"""
Fig 1A–D — Metric explanation panels.

Fig 1A: Raw intensity images (compact vs decompact)
Fig 1B: CV histogram (compact | decompact)
Fig 1C: 1-Gini Lorenz curve (compact | decompact)
Fig 1D: DSI flow panel (normalized image → threshold map → DSI bar)

Optional: Fig 1D alt — DSI 3D surface visualization

Author: Minwoo Kang Ph.D., HR Thiam Lab, Stanford University
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpecFromSubplotSpec as GSFS
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import tifffile as tiff

# ============================================================
# SETTINGS
# ============================================================
DATA_DIR  = os.path.join(os.path.dirname(__file__), "data")
SAVE_SVG  = False
SAVE_DIR  = os.path.join(os.path.dirname(__file__), "outputs")

DSI_THRESHOLD = 0.3
PIXEL_SIZE_UM = 0.098
SCALEBAR_UM   = 5.0
SCALEBAR_PX   = SCALEBAR_UM / PIXEL_SIZE_UM

TITLE_SIZE = 18
LABEL_SIZE = 16
TICK_SIZE  = 14

COLOR_COMPACT   = '#4A90D9'
COLOR_DECOMPACT = '#E5A020'

# ---- Publication style ----
plt.rcParams['font.family']      = 'Arial'
plt.rcParams['axes.linewidth']   = 2.0
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['xtick.direction']  = 'out'
plt.rcParams['ytick.direction']  = 'out'

# ============================================================
# Load images
# ============================================================
img_c  = tiff.imread(os.path.join(DATA_DIR, "From080124_Dish3_Compact_100x100cropped.tif")).astype(np.float64)
img_d  = tiff.imread(os.path.join(DATA_DIR, "From080124_Dish3_Decompact_100x100cropped.tif")).astype(np.float64)
mask_c = tiff.imread(os.path.join(DATA_DIR, "From080124_Dish3_Compact_100x100cropped_mask.tif")) > 0
mask_d = tiff.imread(os.path.join(DATA_DIR, "From080124_Dish3_Decompact_100x100cropped_mask.tif")) > 0

# ============================================================
# Helpers
# ============================================================
def minmax_normalize(vals):
    """Strict min-max normalization to [0,1]."""
    lo, hi = vals.min(), vals.max()
    if hi - lo == 0:
        return np.zeros_like(vals)
    return (vals - lo) / (hi - lo)


def compute_metrics(img, mask):
    """Return raw pixels, normalized pixels, normalized image, and all metric values."""
    vals_raw  = img[mask]
    vals_norm = minmax_normalize(vals_raw)

    # CV (raw intensity)
    mu    = float(np.mean(vals_raw))
    sigma = float(np.std(vals_raw))
    cv    = sigma / mu if mu != 0 else 0.0

    # 1-Gini (normalized intensity)
    x   = np.sort(vals_norm)
    n   = len(x)
    idx = np.arange(1, n + 1)
    s   = float(np.sum(x))
    gini = float(((2 * np.sum(idx * x)) / (n * s)) - ((n + 1) / n)) if s != 0 else 0.0
    one_minus_gini = 1 - gini

    # DSI (normalized intensity)
    dsi = float(np.mean(vals_norm > DSI_THRESHOLD))

    # Normalized image (for DSI map)
    norm_img = np.zeros_like(img, dtype=np.float64)
    norm_img[mask] = vals_norm

    return {
        'vals_raw' : vals_raw,
        'vals_norm': vals_norm,
        'norm_img' : norm_img,
        'mu'       : mu,
        'sigma'    : sigma,
        'cv'       : cv,
        'one_minus_gini': one_minus_gini,
        'dsi'      : dsi,
    }


def add_scalebar(ax, img_shape, scalebar_px, color='white', lw=6, margin_frac=0.05):
    """Draw a scale bar in the lower-right corner."""
    h, w    = img_shape[:2]
    mx, my  = w * margin_frac, h * margin_frac
    x_end   = w - mx
    x_start = x_end - scalebar_px
    y_bar   = h - my
    ax.plot([x_start, x_end], [y_bar, y_bar],
            color=color, linewidth=lw, solid_capstyle='butt')


# Pre-compute metrics
metrics_c = compute_metrics(img_c, mask_c)
metrics_d = compute_metrics(img_d, mask_d)

LABELS = [
    "Early chromatin reorg.\n(Compact)",
    "Late chromatin reorg.\n(Decompact)",
]
DATA   = [metrics_c, metrics_d]
MASKS  = [mask_c, mask_d]
COLORS = [COLOR_COMPACT, COLOR_DECOMPACT]


# ============================================================
# Fig 1A — Raw intensity images
# ============================================================
def make_fig1A(annotate=True):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5),
                             gridspec_kw={'wspace': 0.45})
    rows = [
        ("Early chromatin reorg.\n(Compact)",   img_c, mask_c),
        ("Late chromatin reorg.\n(Decompact)",  img_d, mask_d),
    ]
    for ax, (label, raw_img, mask) in zip(axes, rows):
        im = ax.imshow(raw_img, cmap='gray')
        ax.contour(mask.astype(float), levels=[0.5],
                   colors='yellow', linewidths=2.2, linestyles='--')
        ax.axis('off')
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        if annotate:
            cbar.set_label("Raw intensity", fontsize=12)
            cbar.ax.tick_params(labelsize=10)
            ax.set_title(label, fontsize=TITLE_SIZE)
        else:
            cbar.set_label("")
        add_scalebar(ax, raw_img.shape, SCALEBAR_PX, color='white')
    plt.tight_layout()
    return fig


# ============================================================
# Fig 1B — CV histogram
# ============================================================
def make_fig_cv(annotate=True):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                             gridspec_kw={'wspace': 0.45})
    for ax, m, label, color in zip(axes, DATA, LABELS, COLORS):
        vals = m['vals_raw']
        mu, sigma, cv = m['mu'], m['sigma'], m['cv']

        ax.hist(vals, bins=45, color=color, edgecolor='none', density=False)
        ax.axvline(mu,         linestyle='--', linewidth=2.0, color='#333333')
        ax.axvline(mu - sigma, linestyle=':',  linewidth=2.0, color='#333333')
        ax.axvline(mu + sigma, linestyle=':',  linewidth=2.0, color='#333333')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=TICK_SIZE)

        if annotate:
            ax.set_title(label, fontsize=TITLE_SIZE)
            ax.set_xlabel("Raw intensity (ROI)", fontsize=LABEL_SIZE)
            ax.set_ylabel("Pixel count",         fontsize=LABEL_SIZE)
            ax.text(0.97, 0.95,
                    f"CV = {cv:.3f}\nμ = {mu:.1f}\nσ = {sigma:.1f}",
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=12,
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7, ec='none'))
    plt.tight_layout()
    return fig


# ============================================================
# Fig 1C — 1-Gini Lorenz curve
# ============================================================
def make_fig_1gini(annotate=True):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4),
                             gridspec_kw={'wspace': 0.45})
    for ax, m, label, color in zip(axes, DATA, LABELS, COLORS):
        vals_norm   = m['vals_norm']
        one_minus_g = m['one_minus_gini']

        x       = np.sort(vals_norm)
        cum_int = np.cumsum(x)
        if cum_int[-1] > 0:
            cum_int = cum_int / cum_int[-1]
        cum_px = np.linspace(0, 1, len(cum_int))

        ax.fill_between(cum_px, cum_px, cum_int, alpha=0.45, color='#AAAAAA', label='A')
        ax.fill_between(cum_px, 0, cum_int, alpha=0.30, color=color, label='B')
        ax.plot([0, 1], [0, 1], 'k--', linewidth=2.0)
        ax.plot(cum_px, cum_int, color=color, linewidth=2.6)

        ax.text(0.35, 0.55, 'A', fontsize=16, fontweight='bold',
                color='#555555', transform=ax.transAxes, ha='center')
        ax.text(0.72, 0.22, 'B', fontsize=16, fontweight='bold',
                color=color, transform=ax.transAxes, ha='center')

        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=TICK_SIZE)

        pct_ticks  = [0, 0.25, 0.5, 0.75, 1.0]
        pct_labels = ['0%', '25%', '50%', '75%', '100%']
        ax.set_xticks(pct_ticks); ax.set_xticklabels(pct_labels)
        ax.set_yticks(pct_ticks); ax.set_yticklabels(pct_labels)

        if annotate:
            ax.set_title(label, fontsize=TITLE_SIZE)
            ax.set_xlabel("Cumulative % of pixels (ranked low to high)", fontsize=LABEL_SIZE)
            ax.set_ylabel("Cumulative % of total normalized intensity", fontsize=LABEL_SIZE)
            ax.text(0.50, 0.10, f"1-Gini = {one_minus_g:.3f}",
                    transform=ax.transAxes, ha='center', fontsize=12, fontweight='bold')
    plt.tight_layout()
    return fig


# ============================================================
# Fig 1D — DSI flow panel
# ============================================================
NOT_COUNTED_GRAY   = 0.25
COUNTED_BASE_GRAY  = 0.85
ROI_EDGE_COLOR     = 'yellow'


def make_threshold_class_map(mask, norm_img):
    class_map = np.ones(mask.shape, dtype=float)
    class_map[mask] = NOT_COUNTED_GRAY
    class_map[(norm_img > DSI_THRESHOLD) & mask] = COUNTED_BASE_GRAY
    return class_map


def draw_fraction_bar(ax, dsi_value, color, annotate=True):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    x0, y0, w, h = 0.06, 0.24, 0.88, 0.42
    ax.add_patch(patches.Rectangle(
        (x0, y0), w, h, facecolor=str(NOT_COUNTED_GRAY), edgecolor='black', linewidth=1.6))
    ax.add_patch(patches.Rectangle(
        (x0, y0), w * dsi_value, h, facecolor=color, edgecolor='none'))
    if annotate:
        ax.text(0.50, 0.90, f"DSI = {dsi_value:.3f}",
                ha='center', va='center', fontsize=12, fontweight='bold')
        ax.text(x0,     y0 - 0.08, "0", ha='center', va='top', fontsize=10)
        ax.text(x0 + w, y0 - 0.08, "1", ha='center', va='top', fontsize=10)
        text_x = x0 + w * dsi_value / 2 if dsi_value > 0.10 else x0 + w * dsi_value + 0.04
        text_color = 'white' if dsi_value > 0.18 else 'black'
        ax.text(text_x, y0 + h / 2, f"{dsi_value*100:.1f}%",
                ha='center', va='center', fontsize=10, color=text_color, fontweight='bold')
    ax.axis('off')


def draw_inline_legend(ax, color, annotate=True):
    ax.axis('off')
    if not annotate:
        return
    ax.add_patch(patches.Rectangle(
        (0.03, 0.50), 0.09, 0.25, facecolor=color, edgecolor='none', transform=ax.transAxes))
    ax.text(0.15, 0.62, "Counted",
            transform=ax.transAxes, ha='left', va='center', fontsize=9)
    ax.add_patch(patches.Rectangle(
        (0.50, 0.50), 0.09, 0.25,
        facecolor=str(NOT_COUNTED_GRAY), edgecolor='none', transform=ax.transAxes))
    ax.text(0.62, 0.62, "Not counted",
            transform=ax.transAxes, ha='left', va='center', fontsize=9)


def make_fig1D(annotate=True):
    """Fig 1D: DSI flow panel — normalized image → threshold map → DSI bar."""
    dsi_panels = [
        ("Early chromatin reorg.\n(Compact)",  DATA[0]['norm_img'], mask_c, DATA[0]['dsi'], COLOR_COMPACT),
        ("Late chromatin reorg.\n(Decompact)", DATA[1]['norm_img'], mask_d, DATA[1]['dsi'], COLOR_DECOMPACT),
    ]
    fig = plt.figure(figsize=(13.2, 4.3))
    outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.22)

    for j, (label, norm_img, mask, dsi_val, color) in enumerate(dsi_panels):
        inner = GSFS(2, 3, subplot_spec=outer[0, j],
                     height_ratios=[1.0, 0.24],
                     width_ratios=[1.0, 0.16, 1.0],
                     hspace=0.03, wspace=0.05)

        # Scaled image
        ax_img = fig.add_subplot(inner[0, 0])
        disp = np.zeros_like(norm_img)
        disp[mask] = norm_img[mask]
        im = ax_img.imshow(disp, cmap='gray', vmin=0, vmax=1)
        ax_img.contour(mask.astype(float), levels=[0.5],
                       colors=ROI_EDGE_COLOR, linewidths=1.8, linestyles='--')
        add_scalebar(ax_img, norm_img.shape, SCALEBAR_PX, color='white', lw=6)
        ax_img.axis('off')
        if annotate:
            ax_img.set_title(label, fontsize=14, pad=6)

        # Colorbar with tau line
        cbar = fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.02)
        cbar.set_ticks([0, DSI_THRESHOLD, 1.0])
        if annotate:
            cbar.set_ticklabels(['0', r'$\tau$', '1'])
            cbar.ax.tick_params(labelsize=9)
            cbar.ax.get_yticklabels()[1].set_color('red')
            cbar.ax.get_yticklabels()[1].set_fontweight('bold')
        else:
            cbar.set_ticklabels(['', '', ''])
            cbar.ax.tick_params(length=0)
        cbar.ax.axhline(DSI_THRESHOLD, color='red', linewidth=3.5,
                        xmin=0, xmax=1, clip_on=False, zorder=10)

        # Arrow
        ax_mid = fig.add_subplot(inner[0, 1])
        ax_mid.axis('off')
        if annotate:
            ax_mid.text(0.5, 0.62, "→", ha='center', va='center',
                        fontsize=28, fontweight='bold')
            ax_mid.text(0.5, 0.40,
                        r"$I_{\mathrm{norm}} > \tau$" + "\n" + rf"$\tau$ = {DSI_THRESHOLD}",
                        ha='center', va='center', fontsize=11)

        # Threshold class map
        ax_map = fig.add_subplot(inner[0, 2])
        class_map = make_threshold_class_map(mask, norm_img)
        ax_map.imshow(class_map, cmap='gray', vmin=0, vmax=1)
        overlay = np.zeros((*norm_img.shape, 4), dtype=float)
        r, g, b = plt.matplotlib.colors.to_rgb(color)
        above = (norm_img > DSI_THRESHOLD) & mask
        overlay[above] = [r, g, b, 0.98]
        ax_map.imshow(overlay)
        ax_map.contour(mask.astype(float), levels=[0.5],
                       colors=ROI_EDGE_COLOR, linewidths=1.8, linestyles='--')
        ax_map.axis('off')

        # Bottom blanks
        fig.add_subplot(inner[1, 0]).axis('off')
        fig.add_subplot(inner[1, 1]).axis('off')

        # Legend + DSI bar
        sub_bot = GSFS(2, 1, subplot_spec=inner[1, 2],
                       height_ratios=[0.48, 1.12], hspace=0.04)
        draw_inline_legend(fig.add_subplot(sub_bot[0, 0]), color, annotate)
        draw_fraction_bar(fig.add_subplot(sub_bot[1, 0]), dsi_val, color, annotate)

    if annotate:
        fig.suptitle(r'DSI = fraction of ROI pixels with $I_\mathrm{norm}$ > τ',
                     fontsize=TITLE_SIZE, y=1.02)

    plt.tight_layout()
    return fig


# ============================================================
# Generate & save
# ============================================================
FIGURES = [
    ("Fig1A",       make_fig1A),
    ("Fig1B_CV",    make_fig_cv),
    ("Fig1C_1Gini", make_fig_1gini),
    ("Fig1D_DSI",   make_fig1D),
]

if __name__ == "__main__":
    os.makedirs(SAVE_DIR, exist_ok=True)

    if SAVE_SVG:
        for name, make_fn in FIGURES:
            fig = make_fn(annotate=False)
            fig.savefig(os.path.join(SAVE_DIR, f"{name}.svg"), bbox_inches='tight')
            print(f"Saved: {name}.svg")
            plt.close(fig)
    else:
        for name, make_fn in FIGURES:
            fig = make_fn(annotate=True)
            plt.show()

    print("Done.")
