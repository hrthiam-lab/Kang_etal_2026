#%% =====================================================================

# Shot-noise robustness analysis of CV / 1-Gini / DSI
# (Response to Reviewer 2 / Longo et al.)
#
# All analyses in one script:
#   PART 1 - Nuclear images + added-noise level (compact AND decompact)
#   PART 2 - Single-nucleus metric robustness (% deviation under noise)
#   PART 3 - Compact-vs-decompact separation under noise
#   PART 4 - Per-state metric deviation + direction vs the
#            compact-decompact difference (Longo confounding test)
#
# Same nucleus (cell6), two snapshots: t65 and t159 (compact -> decompact).
# N=500 Poisson shot-noise simulations per condition (mean +/- SD).
# Models photon shot noise only (Poisson); sCMOS read noise NOT added
# (consistent with Longo et al.).

Author: Minwoo Kang Ph.D., HR Thiam Lab, Stanford University
# =====================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from skimage.io import imread

# ---- Output options ----
SAVE_FIG = False          # default: do NOT save (display only)
FIG_FORMAT = "svg"        # vector format when saving
FIG_DPI = 300             # display/raster dpi

# ---- Publication Style ----
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = FIG_DPI
plt.rcParams['savefig.dpi'] = FIG_DPI
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# ---- Parameters ----
DSI_TAU = 0.3
N_REPEATS = 500
# Paper color scheme
COLOR_COMPACT = '#4A90D9'   # blue
COLOR_DECOMPACT = '#E5A020'  # orange
CONV_GAIN = 0.23          # e-/ADU, Kinetix Dynamic Range mode
# Camera bias is included in the local background median (subtracted directly).
rng = np.random.default_rng(42)

# ---- Paths ----
path_dir = r"G:\My Drive\Postdoc_Stanford\Manuscripts\In progress\DSI paper\052726_Nuc_rep_image_from_Nikon_for_Revision"

# Two snapshots of the SAME nucleus (compact -> decompact)
state_files = {
    "t65":  ("3_cell6_spydna_z3_t65_crop.tif",  "3_cell6_spydna_z3_t65_crop_mask.tif"),
    "t159": ("3_cell6_spydna_z3_t159.tif",      "3_cell6_spydna_z3_t159_mask.tif"),
}

imgs, masks = {}, {}
for key, (fimg, fmask) in state_files.items():
    imgs[key] = imread(f"{path_dir}\\{fimg}").astype(np.float64)
    masks[key] = imread(f"{path_dir}\\{fmask}") > 0
    print(f"{key}: img {imgs[key].shape}, mask {masks[key].sum()} px")


def savefig(fig, name):
    """Save figure as vector format only if SAVE_FIG is True."""
    if SAVE_FIG:
        out = f"{path_dir}\\{name}.{FIG_FORMAT}"
        fig.savefig(out, format=FIG_FORMAT, bbox_inches='tight', facecolor='white')
        print(f"Saved: {out}")


# =====================================================================
# Core functions
# =====================================================================
def get_metrics(img, mask, tau=DSI_TAU):
    """Compute CV, 1-Gini, and DSI from masked pixel intensities."""
    vals = img[mask].astype(np.float64)
    if vals.size == 0:
        return 0.0, 0.0, 0.0
    # CV from raw intensities
    cv = np.std(vals) / np.mean(vals) if np.mean(vals) > 0 else 0.0
    # Min-max normalization
    rng_v = vals.max() - vals.min()
    if rng_v == 0:
        return cv, 1.0, 0.0
    vn = (vals - vals.min()) / rng_v
    # DSI
    dsi = np.sum(vn > tau) / len(vn)
    # 1-Gini
    vs = np.sort(vn)
    n = len(vs)
    idx = np.arange(1, n + 1)
    gini = ((2 * np.sum(idx * vs)) / (n * np.sum(vs))) - ((n + 1) / n)
    return cv, 1.0 - gini, dsi


def simulate_photon_noise(img, photon_scale):
    """Poisson shot-noise simulation; photon_scale = max-pixel photon budget.
    The brightest pixel is assigned `photon_scale` expected photons; all other
    pixels scale proportionally. Models photon shot noise only (Poisson)."""
    img_min = img.min()
    img_range = img.max() - img.min()
    if img_range == 0:
        return img.copy()
    img_norm = (img - img_min) / img_range
    expected = np.clip(img_norm * photon_scale, 0.01, None)
    noisy = rng.poisson(expected).astype(np.float64)
    return (noisy / photon_scale) * img_range + img_min


def estimate_max_e(img, mask):
    """Conservative detected-electron estimate for the brightest nuclear pixel."""
    bg_med = np.median(img[~mask])
    return (np.max(img[mask]) - bg_med) * CONV_GAIN


# =====================================================================
# Identify compact vs decompact from original DSI
# =====================================================================
m65, m159 = get_metrics(imgs["t65"], masks["t65"]), get_metrics(imgs["t159"], masks["t159"])
compact_key, decompact_key = ("t65", "t159") if m65[2] < m159[2] else ("t159", "t65")
state_keys = {"compact": compact_key, "decompact": decompact_key}
print(f"\nAssigned: compact = {compact_key}, decompact = {decompact_key} "
      f"(higher DSI = more decompact)")

est_max_e = {st: estimate_max_e(imgs[k], masks[k]) for st, k in state_keys.items()}

# Mean & median detected electrons over NUCLEAR pixels (Longo "average" comparison)
for st, k in state_keys.items():
    bg = np.median(imgs[k][~masks[k]])
    e = (imgs[k][masks[k]] - bg) * CONV_GAIN   # per-pixel detected electrons
    print(f"{st}: mean {e.mean():.0f}, median {np.median(e):.0f}, max {e.max():.0f} e-/px")

# =====================================================================
# Run all simulations once: store per-state, per-metric, per-condition
# =====================================================================
photon_scales = [None, 200, 40]
cond_labels = ["Original\n(no added noise)", "+ Shot noise\n(~200 max photons/px)", "++ Shot noise\n(~40 max photons/px)"]
metric_names = ["CV", "1-Gini", "DSI"]
metric_colors = ['#e74c3c', '#2980b9', '#27ae60']

# Display images: store first realization for each state x condition
display_images = {st: [] for st in state_keys}
# Noise SD accumulator: [state][cond, repeat]
noise_sd = {st: np.zeros((len(photon_scales), N_REPEATS)) for st in state_keys}
# Metric storage: res[state][metric_idx] = [cond, repeat]
res = {st: {mi: np.zeros((len(photon_scales), N_REPEATS)) for mi in range(3)}
       for st in state_keys}

print(f"\nRunning {N_REPEATS} simulations x 2 states x {len(photon_scales)} conditions...")
for st, key in state_keys.items():
    base_img = imgs[key]
    for ci, ps in enumerate(photon_scales):
        for rep in range(N_REPEATS):
            img_sim = base_img.copy() if ps is None else simulate_photon_noise(base_img, ps)
            v = get_metrics(img_sim, masks[key])
            for mi in range(3):
                res[st][mi][ci, rep] = v[mi]
            noise_sd[st][ci, rep] = np.std((img_sim - base_img)[masks[key]])
            if rep == 0:
                display_images[st].append(img_sim)

# Convenience: means/SDs
def msd(st, mi):
    return res[st][mi].mean(axis=1), res[st][mi].std(axis=1)


# =====================================================================
# Console report: PART 2 (robustness), PART 3 (separation), PART 4 (direction)
# =====================================================================
print("\n" + "=" * 80)
print("PART 2 - SINGLE-STATE METRIC ROBUSTNESS (% deviation, orig -> 40-budget)")
print("=" * 80)
for st in state_keys:
    print(f"\n[{st}]  (~{est_max_e[st]:.0f} max e-/px original)")
    for mi, mname in enumerate(metric_names):
        mean, _ = msd(st, mi)
        pct = abs(mean[2] - mean[0]) / mean[0] * 100 if mean[0] != 0 else np.nan
        print(f"  {mname:<7}: {mean[0]:.4f} -> {mean[2]:.4f}  ({pct:.1f}% deviation)")

print("\n" + "=" * 80)
print("PART 3 - COMPACT vs DECOMPACT SEPARATION UNDER NOISE")
print("=" * 80)
for mi, mname in enumerate(metric_names):
    print(f"\n--- {mname} ---")
    cm, _ = msd("compact", mi)
    dm, _ = msd("decompact", mi)
    diff0 = abs(dm[0] - cm[0])
    for ci, clabel in enumerate(cond_labels):
        c, d = res["compact"][mi][ci], res["decompact"][mi][ci]
        pooled = np.sqrt((c.std()**2 + d.std()**2) / 2)
        cohen = abs(d.mean() - c.mean()) / pooled if pooled > 0 else np.inf
        gap = (d.min() - c.max()) if d.mean() > c.mean() else (c.min() - d.max())
        ov = "NO overlap" if gap > 0 else "overlap"
        print(f"  {clabel.replace(chr(10),' '):<22} | compact {c.mean():.3f}+/-{c.std():.3f}"
              f" | decompact {d.mean():.3f}+/-{d.std():.3f} | d={cohen:5.1f} | {ov}")
    worst = max(res["compact"][mi][2].std(), res["decompact"][mi][2].std())
    print(f"  compact-decompact difference (orig) = {diff0:.3f}; "
          f"worst-case noise SD (40) = {worst:.3f}; ratio = {diff0/worst:.1f}")
    
print("\n" + "=" * 80)
print("PART 4 - PER-STATE DEVIATION vs COMPACT-DECOMPACT SEPARATION")
print("=" * 80)
for mi, mname in enumerate(metric_names):
    print(f"\n--- {mname} ---")
    cm, _ = msd("compact", mi)
    dm, _ = msd("decompact", mi)

    # Which state is higher in the noise-free (original) data, and by how much
    diff_signed = (dm[0] - cm[0]) / cm[0] * 100
    diff_dir = "decompact > compact" if dm[0] > cm[0] else "compact > decompact"
    print(f"  Original: compact {cm[0]:.3f} | decompact {dm[0]:.3f} "
          f"(decompact {diff_signed:+.1f}% vs compact; {diff_dir})")

    # Per-state noise-induced deviation (context only)
    for st in ["compact", "decompact"]:
        mean, _ = msd(st, mi)
        dev200 = (mean[1] - mean[0]) / mean[0] * 100
        dev40 = (mean[2] - mean[0]) / mean[0] * 100
        ndir = "increase" if dev40 > 0 else "decrease"
        print(f"  {st:<10}: 200-budget {dev200:+.1f}% | 40-budget {dev40:+.1f}%  "
              f"(noise -> {ndir})")

    # --- Separation (gap) test: the quantity we actually interpret ---
    # Compare the |decompact - compact| gap at 40-budget vs the original gap,
    # rather than inferring direction from one state's movement alone.
    gap_orig = abs(dm[0] - cm[0])
    gap_40 = abs(dm[2] - cm[2])
    gap_change = (gap_40 - gap_orig) / gap_orig * 100 if gap_orig > 0 else np.nan
    if gap_change < -10:
        sep = "REDUCED (conservative: noise shrinks the separation)"
    elif gap_change > 10:
        sep = "INFLATED (caution: noise widens the separation)"
    else:
        sep = "LARGELY PRESERVED"
    print(f"  >> Separation |decompact - compact|: "
          f"orig {gap_orig:.3f} -> 40-budget {gap_40:.3f} "
          f"({gap_change:+.1f}%) -> {sep}")

    # State ranking (which state is higher) preserved at the extreme budget?
    preserved = np.sign(dm[2] - cm[2]) == np.sign(dm[0] - cm[0])
    print(f"  >> Ranking at 40-budget: compact {cm[2]:.3f} vs decompact {dm[2]:.3f} "
          f"-> {'PRESERVED' if preserved else 'REVERSED'}")


# =====================================================================
# FIGURE 1 - Nuclear images (compact & decompact) + added noise level
# =====================================================================
TS, LS, TK, AS = 16, 14, 11, 10

fig1 = plt.figure(figsize=(14, 12))
gs1 = GridSpec(3, 3, figure=fig1, height_ratios=[1, 1, 0.9], hspace=0.45, wspace=0.15)

for row, st in enumerate(["compact", "decompact"]):
    vmin = np.percentile(imgs[state_keys[st]], 1)
    vmax = np.percentile(imgs[state_keys[st]][masks[state_keys[st]]], 99)
    titles = [f"{st.capitalize()} Original\n(~{est_max_e[st]:.0f} max $e^-$/px)",
              f"{st.capitalize()} + Shot noise\n(~200 max photons/px)",
              f"{st.capitalize()} ++ Shot noise\n(~40 max photons/px)"]
    for ci in range(3):
        ax = fig1.add_subplot(gs1[row, ci])
        ax.imshow(display_images[st][ci], cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title(titles[ci], fontsize=TK, fontweight='bold', pad=6)
        ax.axis('off')

# Bottom row: added noise level (mean +/- SD across N=500), grouped by state
ax_noise = fig1.add_subplot(gs1[2, :])
x = np.arange(3)
w = 0.38
for off, st, col in [(-w/2, "compact", COLOR_COMPACT), (w/2, "decompact", COLOR_DECOMPACT)]:
    nm = noise_sd[st].mean(axis=1)
    ns = noise_sd[st].std(axis=1)
    ax_noise.bar(x + off, nm, w, yerr=ns, label=st.capitalize(),
                 color=col, edgecolor='black', linewidth=1.5,
                 capsize=4, error_kw={'linewidth': 1.2})
ax_noise.set_title("Added Noise Level (SD of pixel difference from original)",
                   fontsize=TS - 2, fontweight='bold', pad=8)
ax_noise.set_ylabel("Noise SD (a.u.)", fontsize=LS)
ax_noise.set_xticks(x)
ax_noise.set_xticklabels([l.replace('\n', ' ') for l in cond_labels], fontsize=TK)
ax_noise.tick_params(axis='y', labelsize=TK)
ax_noise.spines['top'].set_visible(False)
ax_noise.spines['right'].set_visible(False)
ax_noise.legend(fontsize=TK, frameon=False)

savefig(fig1, "shot_noise_PART1_images")
plt.show()


# =====================================================================
# FIGURE 2 - Single-state metric robustness (compact & decompact, mean+/-SD)
# (no titles, no in-plot deviation box)
# =====================================================================
# One separate figure per state to avoid x-label / title overlap
for st in ["compact", "decompact"]:
    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
    for mi, (mname, col) in enumerate(zip(metric_names, metric_colors)):
        ax = axes2[mi]
        mean, sd = msd(st, mi)
        bars = ax.bar(x, mean, yerr=sd, width=0.5, color=col, edgecolor='black',
                      linewidth=1.5, alpha=0.85, capsize=5, error_kw={'linewidth': 1.5})
        ax.set_ylabel(f"{mname}  ({st})", fontsize=LS - 1)
        ax.set_xticks(x)
        ax.set_xticklabels([l.replace('\n', ' ') for l in cond_labels],
                           fontsize=TK - 1, rotation=45, ha='right')
        ax.tick_params(axis='y', labelsize=TK)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for bar, mv, sv in zip(bars, mean, sd):
            ax.text(bar.get_x() + bar.get_width()/2, mv + sv + 0.008,
                    f'{mv:.3f}', ha='center', va='bottom', fontsize=AS, fontweight='bold')
        ax.set_ylim(0, max(mean + sd) * 1.30)
    plt.tight_layout()
    savefig(fig2, f"shot_noise_PART2_robustness_{st}")
    plt.show()


# =====================================================================
# FIGURE 3 - Compact vs decompact separation under noise (grouped bars)
# =====================================================================
fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
state_colors = {"compact": COLOR_COMPACT, "decompact": COLOR_DECOMPACT}
for mi, (mname, ax) in enumerate(zip(metric_names, axes3)):
    for off, st in [(-w/2, "compact"), (w/2, "decompact")]:
        mean, sd = msd(st, mi)
        ax.bar(x + off, mean, w, yerr=sd, label=st.capitalize(),
               color=state_colors[st], edgecolor='black', linewidth=1.5,
               capsize=4, error_kw={'linewidth': 1.2})
    ax.set_ylabel(mname, fontsize=LS)
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace('\n', ' ') for l in cond_labels], fontsize=TK - 1, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=TK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # legend removed per request
savefig(fig3, "shot_noise_PART3_separation")
plt.show()


# =====================================================================
# FIGURE 4 - Per-state signed % deviation under noise (direction vs difference)
# =====================================================================
fig4, axes4 = plt.subplots(1, 3, figsize=(15, 5))
for mi, (mname, ax) in enumerate(zip(metric_names, axes4)):
    for st in ["compact", "decompact"]:
        base = res[st][mi][0].mean()
        dev_per_rep = (res[st][mi] - base) / base * 100   # [cond, repeats]
        means = dev_per_rep.mean(axis=1)
        sds = dev_per_rep.std(axis=1)
        ax.errorbar(x, means, yerr=sds, marker='o', markersize=7, linewidth=2,
                    capsize=4, color=state_colors[st], label=st.capitalize())
    ax.axhline(0, color='black', linewidth=1, linestyle='--', alpha=0.5)
    ax.set_ylabel(f"{mname} deviation from original (%)", fontsize=LS - 2)
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace('\n', ' ') for l in cond_labels], fontsize=TK - 1, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=TK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if mi == 2:
        ax.legend(fontsize=TK, frameon=False, loc='best')
savefig(fig4, "shot_noise_PART4_deviation")
plt.show()

print("\nDone. Set SAVE_FIG=True to export vector (SVG) figures.")
# Output folders mirror the live layout, created next to each cell TIFF:
#   TIFF/1.ome/1.ome_cell_1_SPYDNA/        (event signal frames)
#   TIFF/1.ome/1.ome_cell_1_SPYDNA_mask/   (masks)