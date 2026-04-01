"""
Fig 2: Live-cell NETosis tracking — CV, 1-Gini, DSI
=====================================================
Row 1 (B–D): Time-series (mean ± SD, normalized time) + p-value on twinx
Row 2 (B–D): Trajectory-level summary (per-cell mean) + Mann-Whitney U

Metric order: CV | 1-Gini | DSI
Statistical tests:
  - Time-series: Brunner-Munzel + BH FDR correction
  - Trajectory-level: Mann-Whitney U (two-sided)

Author: Minwoo Kang Ph.D., Hawa Thiam Lab, Stanford University
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import brunnermunzel, mannwhitneyu, linregress
from statsmodels.stats.multitest import multipletests

# ============================================================
# SETTINGS
# ============================================================
DATA_DIR  = os.path.join(os.path.dirname(__file__), "data")
SAVE_SVG  = False
SAVE_DIR  = os.path.join(os.path.dirname(__file__), "outputs")

DATA_PATH = os.path.join(DATA_DIR, "NETing_NonNETing_SPYDNA_features.csv")

# Metric order: CV, 1-Gini, DSI
METRICS = [
    ('intensity_std / intensity_mean', 'CV',     'CV (σ/μ)'),
    ('SDI',                            '1-Gini', '1-Gini'),
    ('DSI',                            'DSI',    'DSI'),
]

COLOR_NETING     = '#E55C5C'
COLOR_NON_NETING = '#4A90D9'
COLOR_ALPHA      = 0.25

COMMON_TIMES = np.linspace(0, 1, 241)
A_MIN, A_MAX = 21.71, 131.62

TITLE_SIZE = 22
LABEL_SIZE = 20
TICK_SIZE  = 22

# ---- Publication style ----
rcParams['font.family']      = 'Arial'
rcParams['axes.linewidth']   = 2.0
rcParams['xtick.major.width'] = 2.0
rcParams['ytick.major.width'] = 2.0
rcParams['xtick.direction']  = 'out'
rcParams['ytick.direction']  = 'out'

# ============================================================
# Load & prep
# ============================================================
data = pd.read_csv(DATA_PATH)

if 'area' in data.columns:
    inv = (1.0 / data['area'] - 1.0 / A_MAX) / (1.0 / A_MIN - 1.0 / A_MAX)
    data['inv_area_norm'] = np.clip(inv, 0.0, 1.0)

data['unique_cell_id'] = (data['image_stack'].astype(str) + '_' +
                          data['cell_label'].astype(str))


def normalize_time(group):
    tmin, tmax = group.min(), group.max()
    return (group - tmin) / (tmax - tmin) if tmax > tmin else group * 0.0


data['normalized_time'] = data.groupby('unique_cell_id')['time_point'].transform(normalize_time)

n_neting     = data.loc[data['cell_status'] == 'NETing',     'unique_cell_id'].nunique()
n_non_neting = data.loc[data['cell_status'] == 'non_NETing', 'unique_cell_id'].nunique()
print(f"NETing cells    : {n_neting}")
print(f"non-NETing cells: {n_non_neting}")

# ============================================================
# Helpers
# ============================================================
def evaluate_feature(df, expression):
    """Return Series for column name or arithmetic expression."""
    if expression in df.columns:
        return df[expression]
    try:
        return df.eval(expression, engine='python')
    except Exception as e:
        raise ValueError(f"Cannot evaluate '{expression}': {e}")


def interpolate_cell(cell_df, time_grid, col):
    """Linear interpolation of one cell onto common time grid."""
    if cell_df.empty:
        return np.full(time_grid.shape, np.nan)
    cell_df = cell_df.sort_values('normalized_time')
    x = cell_df['normalized_time'].to_numpy()
    y = cell_df[col].to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 2:
        return np.full(time_grid.shape, np.nan)
    return np.interp(time_grid, x, y, left=np.nan, right=np.nan)


def build_matrix(df, time_grid, col, status):
    """Build (n_cells x n_timepoints) matrix for one group."""
    rows = []
    sub  = df[df['cell_status'] == status]
    for uid in sub['unique_cell_id'].unique():
        rows.append(interpolate_cell(sub[sub['unique_cell_id'] == uid], time_grid, col))
    return np.array(rows, dtype=float)


def bm_fdr(mat_a, mat_b, alpha=0.05, min_n=3):
    """Brunner-Munzel per time point + BH FDR correction."""
    T     = mat_a.shape[1]
    pvals = np.full(T, np.nan)
    for i in range(T):
        x = mat_a[:, i][np.isfinite(mat_a[:, i])]
        y = mat_b[:, i][np.isfinite(mat_b[:, i])]
        if len(x) < min_n or len(y) < min_n:
            continue
        try:
            _, p = brunnermunzel(x, y, alternative='two-sided', distribution='t')
            pvals[i] = p
        except Exception:
            pass
    p_adj = np.full(T, np.nan)
    ok = np.isfinite(pvals)
    if ok.any():
        _, padj_ok, _, _ = multipletests(pvals[ok], method='fdr_bh', alpha=alpha)
        p_adj[ok] = padj_ok
    return p_adj


def per_cell_summary(df, time_grid, col, status):
    """Summarize each cell trajectory as scalar (mean over normalized time)."""
    records = []
    sub = df[df['cell_status'] == status]
    for uid in sub['unique_cell_id'].unique():
        vals = interpolate_cell(sub[sub['unique_cell_id'] == uid], time_grid, col)
        ok   = np.isfinite(vals)
        if ok.sum() < 5:
            continue
        records.append({
            'unique_cell_id': uid,
            'status'        : status,
            'traj_mean'     : float(np.mean(vals[ok])),
            'traj_slope'    : float(linregress(time_grid[ok], vals[ok]).slope),
        })
    return pd.DataFrame(records)


# ============================================================
# Pre-compute all metrics
# ============================================================
results = {}

for expr, label, ylabel in METRICS:
    print(f"  Processing: {label} ...")
    data['_feat_raw'] = evaluate_feature(data, expr)
    data['_feat_raw'] = data['_feat_raw'].replace([np.inf, -np.inf], np.nan)

    net_raw = build_matrix(data, COMMON_TIMES, '_feat_raw', 'NETing')
    non_raw = build_matrix(data, COMMON_TIMES, '_feat_raw', 'non_NETing')

    df_net = per_cell_summary(data, COMMON_TIMES, '_feat_raw', 'NETing')
    df_non = per_cell_summary(data, COMMON_TIMES, '_feat_raw', 'non_NETing')

    results[expr] = {
        'label'    : label,
        'ylabel'   : ylabel,
        'net_mean' : np.nanmean(net_raw, axis=0),
        'net_std'  : np.nanstd(net_raw,  axis=0),
        'non_mean' : np.nanmean(non_raw, axis=0),
        'non_std'  : np.nanstd(non_raw,  axis=0),
        'p_raw'    : bm_fdr(net_raw, non_raw),
        'net_traj' : df_net,
        'non_traj' : df_non,
    }

print("Done computing.\n")

# ============================================================
# Row 1: Time-series panel
# ============================================================
def draw_timeseries(ax, expr, annotate=True, show_legend=False, show_pval_label=False):
    r = results[expr]
    t = COMMON_TIMES

    ax.plot(t, r['net_mean'],  color=COLOR_NETING,     linewidth=2.5, label='NETing')
    ax.fill_between(t, r['net_mean'] - r['net_std'], r['net_mean'] + r['net_std'],
                    color=COLOR_NETING, alpha=COLOR_ALPHA)
    ax.plot(t, r['non_mean'],  color=COLOR_NON_NETING, linewidth=2.5, label='non-NETing')
    ax.fill_between(t, r['non_mean'] - r['non_std'], r['non_mean'] + r['non_std'],
                    color=COLOR_NON_NETING, alpha=COLOR_ALPHA)

    ax2 = ax.twinx()
    ax2.plot(t, r['p_raw'], color='#AAAAAA', linewidth=1.0, linestyle='--',
             alpha=0.7, label='p-value (FDR)')
    ax2.axhline(0.05, color='#888888', linewidth=1.0, linestyle='-', label='p=0.05')
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_yticks([0, 0.5, 1.0])
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(labelsize=TICK_SIZE)

    ax.spines['top'].set_visible(False)
    ax.tick_params(labelsize=TICK_SIZE)

    if annotate:
        ax.set_xlabel('Normalized time', fontsize=LABEL_SIZE)
        ax.set_ylabel(r['ylabel'], fontsize=LABEL_SIZE)
        ax.set_title(r['label'], fontsize=TITLE_SIZE)
        if show_pval_label:
            ax2.set_ylabel('p-value (BH-FDR)', fontsize=LABEL_SIZE - 2)
        if show_legend:
            lines1, labs1 = ax.get_legend_handles_labels()
            lines2, labs2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labs1 + labs2,
                      fontsize=12, loc='upper center',
                      bbox_to_anchor=(0.5, 1.30), ncol=2)
    else:
        ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')
        ax2.set_ylabel('')


# ============================================================
# Row 2: Trajectory summary panel
# ============================================================
MEAN_LINE_W = 0.18
CAP_W       = 0.08


def draw_trajectory(ax, expr, annotate=True):
    r        = results[expr]
    net_vals = r['net_traj']['traj_mean'].values
    non_vals = r['non_traj']['traj_mean'].values

    groups = [(1, non_vals, COLOR_NON_NETING, 'non-NETing'),
              (2, net_vals, COLOR_NETING,     'NETing')]

    np.random.seed(42)
    for pos, vals, col, lbl in groups:
        mean = np.mean(vals)
        sd   = np.std(vals)

        jitter = np.random.uniform(-0.08, 0.08, len(vals))
        ax.scatter(pos + jitter, vals, color=col, alpha=0.6, s=45, zorder=2, label=lbl)

        ax.plot([pos, pos], [mean - sd, mean + sd],
                color='black', linewidth=2.5, zorder=3, solid_capstyle='butt')
        ax.plot([pos - CAP_W, pos + CAP_W], [mean + sd, mean + sd],
                color='black', linewidth=2.5, zorder=3)
        ax.plot([pos - CAP_W, pos + CAP_W], [mean - sd, mean - sd],
                color='black', linewidth=2.5, zorder=3)
        ax.plot([pos - MEAN_LINE_W, pos + MEAN_LINE_W], [mean, mean],
                color='black', linewidth=2.5, zorder=4, solid_capstyle='butt')

    _, p = mannwhitneyu(net_vals, non_vals, alternative='two-sided')
    y_min   = min(np.min(net_vals), np.min(non_vals))
    y_max   = max(np.max(net_vals), np.max(non_vals))
    y_range = y_max - y_min
    y_top   = y_max + y_range * 0.15
    y_text  = y_top + y_range * 0.03

    if annotate:
        ax.plot([1, 2], [y_top, y_top], 'k-', linewidth=1.5)
        stars = ('***' if p < 0.001 else '**' if p < 0.01
                 else '*' if p < 0.05 else 'ns')
        ax.text(1.5, y_text, f"{stars}  (p={p:.2e})",
                ha='center', va='bottom', fontsize=12)

    y_bottom = y_min - y_range * 0.1
    ax.set_ylim(y_bottom, y_text + y_range * 0.15)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=TICK_SIZE)

    if annotate:
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['non-NETing', 'NETing'], fontsize=LABEL_SIZE - 1)
        ax.set_ylabel(f"Mean {r['ylabel']}\n(trajectory avg)", fontsize=LABEL_SIZE)
        ax.text(0.98, 0.02, f'n={len(non_vals)} / {len(net_vals)}',
                transform=ax.transAxes, ha='right', fontsize=10, color='gray')
    else:
        ax.set_xticks([1, 2])
        ax.set_xticklabels([])
        ax.set_ylabel('')


# ============================================================
# Combined figure
# ============================================================
def make_fig(annotate=True):
    fig, axes = plt.subplots(2, 3, figsize=(26, 14),
                             gridspec_kw={'hspace': 0.5, 'wspace': 0.45})

    for col, (expr, _, _) in enumerate(METRICS):
        draw_timeseries(axes[0, col], expr,
                        annotate       =annotate,
                        show_legend    =(col == 0),
                        show_pval_label=(col == 2))
        draw_trajectory(axes[1, col], expr, annotate=annotate)

    if annotate:
        fig.text(0.01, 0.01,
                 f'NETing n={n_neting}  |  non-NETing n={n_non_neting}',
                 ha='left', va='bottom', fontsize=12, color='gray')
    plt.tight_layout()
    return fig


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    os.makedirs(SAVE_DIR, exist_ok=True)

    if SAVE_SVG:
        fig = make_fig(annotate=False)
        fig.savefig(os.path.join(SAVE_DIR, "Fig2_combined.svg"), bbox_inches='tight')
        print("Saved: Fig2_combined.svg")
        plt.close(fig)
    else:
        fig = make_fig(annotate=True)
        plt.show()

    print("Done.")
