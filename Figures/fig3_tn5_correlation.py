"""
Fig 3: Tn5 accessibility correlation
======================================
Panel layout (1 row × 3 columns):
  [A] CV vs Tn5          — Spearman ρ, linear fit + inset zoom
  [B] 1-Gini vs Tn5      — Spearman ρ, linear fit + inset zoom
  [C] DSI vs Tn5         — Spearman ρ, linear fit + inset zoom

Statuses shown: Bef Nuclear rounding, Nuclear rounding

Author: Minwoo Kang Ph.D., Hawa Thiam Lab, Stanford University
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from pathlib import Path
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr, pearsonr

# ============================================================
# SETTINGS
# ============================================================
DATA_DIR  = os.path.join(os.path.dirname(__file__), "data")
SAVE_SVG  = False
SAVE_DIR  = os.path.join(os.path.dirname(__file__), "outputs")

X_PATH = os.path.join(DATA_DIR, "pooled_DAPI.csv")
Y_PATH = os.path.join(DATA_DIR, "pooled_Tn5.csv")

SELECTED_STATUSES_CSV = ["Bef Nuclear rounding", "Nuclear rounding"]

STATUS_DISPLAY = {
    "Bef Nuclear rounding": "Bef Nuclear rounding",
    "Nuclear rounding"    : "Nuclear rounding",
}

COLOR_MAP = {
    "Bef Nuclear rounding": "#fcad91",
    "Nuclear rounding"    : "#fb6a4a",
}

CORR_METRICS = [
    ("intensity_std / intensity_mean", "CV"),
    ("SDI",                            "1-Gini"),
    ("DSI",                            "DSI"),
]

INSET_X_PERCENTILE = (5, 95)
INSET_LOC = ['upper right', 'upper left', 'upper left']

LABEL_SIZE  = 22
TICK_SIZE   = 28
SPINE_W     = 1
TICK_W      = 2.5
TICK_LEN    = 5
SCATTER_S   = 50
SCATTER_A   = 0.55

# ---- Publication style ----
rcParams['font.family']  = 'Arial'
rcParams['figure.dpi']   = 300

# ============================================================
# Helper functions
# ============================================================
def standardize_status(series):
    s = series.astype(str).str.strip().str.lower()
    mapper = {
        "nuclear rounding"         : "Nuclear rounding",
        "nuclear  rounding"        : "Nuclear rounding",
        "nuclear-rounding"         : "Nuclear rounding",
        "nuclear rounding (early)" : "Nuclear rounding",
        "before nuclear rounding"  : "Bef Nuclear rounding",
        "bef nuclear rounding"     : "Bef Nuclear rounding",
        "before nuclear rupture"   : "Bef Nuclear rounding",
        "bef nuclear rupture"      : "Bef Nuclear rounding",
        "before nulcear rupture"   : "Bef Nuclear rounding",
        "ne rupture"               : "NE rupture",
        "ne-rupture"               : "NE rupture",
        "ne  rupture"              : "NE rupture",
        "non-stimulated"           : "Non-stimulated",
        "non stimulated"           : "Non-stimulated",
        "netosis/pm rupture"       : "NETosis/PM rupture",
    }
    out = s.map(mapper)
    return series.where(out.isna(), out)


def _norm_path(s):
    s = str(s).replace("\\", "/")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def extract_exp_pos_key(fullpath):
    s = _norm_path(fullpath)
    parts = [p for p in s.split("/") if p]
    for i in range(len(parts) - 1, -1, -1):
        if re.fullmatch(r"Position[_-]?\d+", parts[i], flags=re.IGNORECASE):
            pos = re.sub(r"Position[-_]?(\d+)", r"Position_\1", parts[i], flags=re.IGNORECASE)
            exp = parts[i-1] if i-1 >= 0 else ""
            return f"{exp}/{pos}" if exp else pos
    return None


def extract_trailing_index(fullpath):
    stem = Path(_norm_path(fullpath)).stem
    m = re.search(r"(\d+)$", stem)
    return int(m.group(1)) if m else None


def ensure_integrated_side(df, poi, tag, out_col):
    if out_col in df.columns:
        return
    area_candidates = [f"area_{poi}_{tag}", f"area_{tag}", f"area_{poi}", "area"]
    mean_candidates = [
        f"intensity_mean_{poi}_{tag}", f"intensity_mean_{tag}", f"intensity_mean_{poi}", "intensity_mean",
        f"mean_intensity_{poi}_{tag}", f"mean_intensity_{tag}", f"mean_intensity_{poi}", "mean_intensity",
    ]
    area_col = next((c for c in area_candidates if c in df.columns), None)
    mean_col = next((c for c in mean_candidates if c in df.columns), None)
    if area_col and mean_col:
        with np.errstate(invalid='ignore'):
            df[out_col] = (pd.to_numeric(df[area_col], errors='coerce') *
                           pd.to_numeric(df[mean_col], errors='coerce'))


_ALLOWED_FUNCS = {
    "abs": np.abs, "sqrt": np.sqrt, "log": np.log, "log10": np.log10,
    "exp": np.exp, "clip": np.clip, "where": np.where,
    "min": np.minimum, "max": np.maximum, "pow": np.power,
}
_TOKEN_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _resolve_col(df, base, poi, tag):
    for c in [f"{base}_{poi}_{tag}", f"{base}_{tag}", f"{base}_{poi}", base]:
        if c in df.columns:
            return c
    return None


def build_values(df, feature, poi, tag):
    """Resolve a feature expression to numpy array + display label."""
    f = feature.strip()

    if f.lower() in {"integrated_intensity", "integrated", "int_int"}:
        out_col = f"integrated_intensity_{tag}"
        ensure_integrated_side(df, poi=poi, tag=tag, out_col=out_col)
        for c in [out_col, f"integrated_intensity_{poi}_{tag}",
                  f"integrated_intensity_{tag}", f"integrated_intensity_{poi}",
                  "integrated_intensity"]:
            if c in df.columns:
                return pd.to_numeric(df[c], errors='coerce').to_numpy(), f"integrated_intensity ({poi})"
        raise KeyError(f"Cannot resolve integrated_intensity on side {tag}")

    if not any(op in f for op in "+-*/()"):
        col = _resolve_col(df, f, poi, tag)
        if col is None:
            raise KeyError(f"Column '{f}' not found for side {tag}")
        return pd.to_numeric(df[col], errors='coerce').to_numpy(), f"{f} ({poi})"

    tokens = set(_TOKEN_RE.findall(f))
    local_dict = dict(_ALLOWED_FUNCS)
    for t in tokens:
        if t in _ALLOWED_FUNCS:
            continue
        col = _resolve_col(df, t, poi, tag)
        if col is None:
            raise KeyError(f"Unresolved symbol '{t}' for side {tag}")
        local_dict[t] = pd.to_numeric(df[col], errors='coerce')
    vals = pd.eval(f, local_dict=local_dict, engine='python')
    return pd.to_numeric(vals, errors='coerce').to_numpy(), f"{f} (DAPI)"


# ============================================================
# Load & match data
# ============================================================
X = pd.read_csv(X_PATH)
Y = pd.read_csv(Y_PATH)

for df in (X, Y):
    df['cell_status'] = standardize_status(df['cell_status']) if 'cell_status' in df.columns else np.nan
    df['exp_pos_key'] = df['file_name'].map(extract_exp_pos_key)
    df['num_id']      = df['file_name'].map(extract_trailing_index)

X['rank_in_key'] = X.groupby(['exp_pos_key', 'num_id']).cumcount()
Y['rank_in_key'] = Y.groupby(['exp_pos_key', 'num_id']).cumcount()

merged = pd.merge(X, Y, on=['exp_pos_key', 'num_id', 'rank_in_key'],
                  how='inner', suffixes=('_X', '_Y'))
merged = merged[merged['cell_status_X'].isin(SELECTED_STATUSES_CSV)].copy()

print(f"[Match] {len(merged)} paired cells after status filter")
print(merged['cell_status_X'].value_counts().to_string())

# ============================================================
# Build correlation data
# ============================================================
tn5_vals, _ = build_values(merged, "integrated_intensity", "Tn5", tag="Y")

corr_data = []
for feat_csv, display_label in CORR_METRICS:
    xv, _ = build_values(merged, feat_csv, "DAPI", tag="X")
    mask   = np.isfinite(xv) & np.isfinite(tn5_vals)
    df_sub = merged.loc[mask].copy()
    df_sub["__X"] = xv[mask]
    df_sub["__Y"] = tn5_vals[mask]

    xv_clean = df_sub["__X"].to_numpy()
    yv_clean = df_sub["__Y"].to_numpy()

    model  = LinearRegression().fit(xv_clean.reshape(-1, 1), yv_clean)
    y_pred = model.predict(xv_clean.reshape(-1, 1))
    rho, p_rho = spearmanr(xv_clean, yv_clean)
    r_p, p_p   = pearsonr(xv_clean, yv_clean)

    x_lo = np.percentile(xv_clean, INSET_X_PERCENTILE[0])
    x_hi = np.percentile(xv_clean, INSET_X_PERCENTILE[1])

    corr_data.append({
        'display_label': display_label,
        'df'           : df_sub,
        'xv'           : xv_clean,
        'yv'           : yv_clean,
        'y_pred'       : y_pred,
        'model'        : model,
        'rho'          : rho,
        'p_rho'        : p_rho,
        'inset_xlim'   : (x_lo, x_hi),
    })
    print(f"  {display_label:8s}  n={mask.sum():3d} | "
          f"Spearman ρ={rho:+.3f} (p={p_rho:.2e}) | "
          f"Pearson  r={r_p:+.3f} (p={p_p:.2e})")

# ============================================================
# Figure
# ============================================================
def make_fig(annotate=True):
    fig, axes = plt.subplots(1, 3, figsize=(22, 7),
                             gridspec_kw={'wspace': 0.55})

    for i, cd in enumerate(corr_data):
        ax     = axes[i]
        df_sub = cd['df']
        model  = cd['model']
        rho    = cd['rho']
        p_rho  = cd['p_rho']
        label  = cd['display_label']

        # Main scatter
        for status_csv in SELECTED_STATUSES_CSV:
            sub = df_sub[df_sub['cell_status_X'] == status_csv]
            if sub.empty:
                continue
            disp = STATUS_DISPLAY.get(status_csv, status_csv)
            ax.scatter(sub["__X"], sub["__Y"],
                       s=SCATTER_S, alpha=SCATTER_A,
                       c=COLOR_MAP.get(status_csv, 'gray'),
                       label=f"{disp} (n={len(sub)})")

        ax.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
        ax.yaxis.get_offset_text().set_fontsize(TICK_SIZE)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(SPINE_W)
        ax.spines['bottom'].set_linewidth(SPINE_W)
        ax.tick_params(labelsize=TICK_SIZE, width=TICK_W, length=TICK_LEN)

        main_ylim = ax.get_ylim()

        if annotate:
            ax.set_xlabel(label, fontsize=LABEL_SIZE, labelpad=8)
            ax.set_ylabel("Tn5\nintegrated intensity", fontsize=LABEL_SIZE, labelpad=8)
            ax.set_title(f"Spearman's ρ = {rho:.3f}\n(p={p_rho:.2e})",
                         fontsize=12, pad=10)
            if i == 2:
                ax.legend(fontsize=10, frameon=False, loc='lower right')

        # Inset zoom
        if INSET_LOC[i] == 'upper left':
            axins = inset_axes(ax, width="45%", height="45%",
                               bbox_to_anchor=(0.08, 0.0, 1.0, 1.0),
                               bbox_transform=ax.transAxes,
                               loc='upper left', borderpad=1.5)
        else:
            axins = inset_axes(ax, width="45%", height="45%",
                               loc=INSET_LOC[i], borderpad=1.5)

        axins.patch.set_facecolor('white')
        axins.patch.set_alpha(0.85)

        for status_csv in SELECTED_STATUSES_CSV:
            sub = df_sub[df_sub['cell_status_X'] == status_csv]
            if sub.empty:
                continue
            axins.scatter(sub["__X"], sub["__Y"],
                          s=SCATTER_S * 0.6, alpha=SCATTER_A + 0.1,
                          c=COLOR_MAP.get(status_csv, 'gray'), edgecolors='none')

        inset_xlim = cd['inset_xlim']
        x_fit_inset = np.linspace(inset_xlim[0], inset_xlim[1], 100)
        y_fit_inset = model.predict(x_fit_inset.reshape(-1, 1))
        axins.plot(x_fit_inset, y_fit_inset, color='black', linewidth=1.2,
                   linestyle='-', zorder=5)

        axins.set_xlim(inset_xlim)
        axins.set_ylim(main_ylim)

        axins.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
        axins.yaxis.get_offset_text().set_fontsize(12)
        axins.tick_params(axis='y', labelsize=16, width=1.0, length=3)
        axins.tick_params(axis='x', labelsize=16, width=1.0, length=3)

        axins.spines['top'].set_visible(True)
        axins.spines['right'].set_visible(True)
        for spine in axins.spines.values():
            spine.set_linewidth(1.2)
            spine.set_edgecolor('#666666')

        if INSET_LOC[i] == 'upper left':
            mark_inset(ax, axins, loc1=1, loc2=3,
                       fc="none", ec="#999999", linewidth=1.0, linestyle='--')
        else:
            mark_inset(ax, axins, loc1=2, loc2=4,
                       fc="none", ec="#999999", linewidth=1.0, linestyle='--')

    plt.tight_layout()
    return fig


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    os.makedirs(SAVE_DIR, exist_ok=True)

    if SAVE_SVG:
        fig = make_fig(annotate=False)
        fig.savefig(os.path.join(SAVE_DIR, "Fig3_scatter.svg"), bbox_inches='tight')
        print("Saved: Fig3_scatter.svg")
        plt.close(fig)
    else:
        fig = make_fig(annotate=True)
        plt.show()

    print("Done.")
