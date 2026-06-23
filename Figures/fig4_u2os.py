#%% =====================================================================
"""
Fig 4: U2OS mitotic chromatin compaction validation
====================================================
Panel layout:
  [A] Representative U2OS images
      — Interphase (decompact chromatin)
      — Mitotic entry (compact chromatin)
      — Post-mitotic daughter nuclei (decompact chromatin)
      — Panel A is assembled separately in Adobe Illustrator.

  [B] CV across mitotic stages
  [C] 1-Gini across mitotic stages
  [D] DSI across mitotic stages

Data:
  - Multiple independent U2OS live-imaging experiments are pooled.
  - The same cell is measured at three manually selected in-focus stages.
  - Raw cell_id values may overlap between experiments; therefore,
    pairing uses a globally unique identifier: experiment + cell_id.
  - Chromosome segregation frames are excluded by default because
    chromosome movement and mitotic rounding frequently displace the
    DNA-containing region from the single focal plane.

Statistics:
  - Friedman test for the paired omnibus stage comparison.
  - Kendall's W for the omnibus effect size.
  - Paired Wilcoxon signed-rank tests for all pairwise comparisons.
  - Holm correction across the three pairwise comparisons.
  - Matched-pairs rank-biserial correlation for pairwise effect size.

Visual encoding:
  - Mitotic entry / compact chromatin: blue (#4A90D9).
  - Interphase and post-mitotic / decompact chromatin: dark gray tones.
  - Paired-cell connecting lines: light gray.
  - Point size and summary-bar styling match revised Figure 2.

Output:
  - Generates Figure 4B-D as one 1 × 3 vector figure.
  - SVG text remains editable in Adobe Illustrator.
  - Final physical dimensions are specified in millimeters.
"""
# =====================================================================

import os
import glob
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests


# ==============================================================
# Publication style
# ==============================================================

def mm_to_inch(mm):
    """Convert millimeters to inches for matplotlib figsize."""
    return mm / 25.4


rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
    ],

    # Final manuscript font sizes
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 8,

    # Final manuscript line widths
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "lines.linewidth": 1.0,

    # Keep vector text editable in Illustrator
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ==============================================================
# Settings
# ==============================================================

PARENT_DIR = r"E:\Live\U2OS_mitosis_for_revision_Nulceus"
CSV_GLOB = "*mitotic_metrics.csv"

SAVE_DIR = (
    r"G:\My Drive\Postdoc_Stanford"
    r"\Manuscripts\In progress\DSI paper"
)

# Review the figure first.
SAVE_SVG = False
SAVE_PDF = False
SAVE_STATS_CSV = False

OUTPUT_BASENAME = "Fig4_U2OS_mitotic_metrics"

# Final physical dimensions of panels B-D
FINAL_FIG_WIDTH_MM = 180
FINAL_FIG_HEIGHT_MM = 72

# Segregation is excluded because single-Z imaging is focus-compromised.
DROP_SEGREGATION = True

SHOW_PAIRED_LINES = True

# Keep statistics in the caption rather than the panel titles by default.
SHOW_STATS_IN_TITLE = False

METRICS = [
    ("CV", "CV"),
    ("1-Gini", "1-Gini"),
    ("DSI", "DSI"),
]

EVENT_MAP = {
    "Interphase (pre-mitotic)": "Interphase",
    "Mitotic condensation (prometaphase–metaphase)": "Mitotic entry",
    "Chromosome segregation (anaphase–telophase)": "Segregation",
    "Post-mitotic (daughter nuclei)": "Post-mitotic",
}

DISPLAY_LABEL = {
    "Interphase": "Interphase\n(decompact)",
    "Mitotic entry": "Mitotic entry\n(compact)",
    "Segregation": "Segregation",
    "Post-mitotic": "Post-mitotic\n(decompact)",
}

# Compact chromatin is blue.
# Decompact states are dark gray tones.
COLOR = {
    "Interphase": "#555555",
    "Mitotic entry": "#4A90D9",
    "Segregation": "#8A6FA8",
    "Post-mitotic": "#777777",
}

# Paired lines are lighter than all stage markers.
COLOR_PAIRED_LINE = "#D0D0D0"

# Final font sizes in points
TITLE_SIZE = 10
AXIS_LABEL_SIZE = 10
X_TICK_LABEL_SIZE = 10
Y_TICK_LABEL_SIZE = 8
SIGNIFICANCE_FONT_SIZE = 8

# Final line widths in points
AXIS_LINE_WIDTH = 1.0
PAIRED_LINE_WIDTH = 0.35
SD_LINE_WIDTH = 1.0
MEAN_LINE_WIDTH = 1.0
BRACKET_LINE_WIDTH = 0.8

# Scatter styling matched to revised Figure 2
SCATTER_SIZE = 10
SCATTER_ALPHA = 0.60
PAIRED_LINE_ALPHA = 0.75

# Mean ± SD geometry in x-axis data coordinates
MEAN_LINE_HALF_WIDTH = 0.22
SD_CAP_HALF_WIDTH = 0.09

# Jitter
JITTER_WIDTH = 0.09
JITTER_SEED = 42


# ==============================================================
# Load and combine CSV files
# ==============================================================

def extract_experiment_token(path):
    """Extract the six-digit experiment identifier from a path."""
    for path_component in re.split(r"[\\/]+", path):
        match = re.match(r"(\d{6})", path_component)

        if match:
            return match.group(1)

    return os.path.basename(
        os.path.dirname(path)
    )


csv_paths = sorted(
    glob.glob(
        os.path.join(
            PARENT_DIR,
            "**",
            CSV_GLOB,
        ),
        recursive=True,
    )
)

if not csv_paths:
    raise FileNotFoundError(
        f"No files matching '{CSV_GLOB}' were found under:\n"
        f"{PARENT_DIR}"
    )


frames = []

print("Loading metric CSVs:")

for csv_path in csv_paths:
    experiment_df = pd.read_csv(csv_path)

    experiment_df["experiment"] = extract_experiment_token(
        csv_path
    )

    frames.append(experiment_df)

    print(
        f"  {len(experiment_df):>4d} rows  "
        f"exp={experiment_df['experiment'].iloc[0]:<8}  "
        f"<- {csv_path}"
    )


df = pd.concat(
    frames,
    ignore_index=True,
)

df["ev"] = df["event"].map(EVENT_MAP)

unmapped_events = (
    df.loc[df["ev"].isna(), "event"]
    .dropna()
    .unique()
)

if len(unmapped_events) > 0:
    print("\n[WARN] The following event labels were not mapped:")

    for event_name in unmapped_events:
        print(f"  - {event_name}")


# Construct a globally unique cell identifier.
df["uid"] = (
    df["experiment"].astype(str)
    + "_"
    + df["cell_id"].astype(str)
)

n_experiments = df["experiment"].nunique()
raw_cell_count = df["cell_id"].nunique()
unique_cell_count = df["uid"].nunique()

collision_text = (
    "COLLISION present"
    if unique_cell_count != raw_cell_count
    else "no collision"
)

print(f"\n{n_experiments} experiment(s) combined.")

print(
    f"Unique cells by uid: {unique_cell_count} "
    f"(raw cell_id would have given {raw_cell_count} "
    f"-> {collision_text})"
)


order = [
    "Interphase",
    "Mitotic entry",
    "Segregation",
    "Post-mitotic",
]

if DROP_SEGREGATION:
    order = [
        event
        for event in order
        if event != "Segregation"
    ]

df = df[df["ev"].isin(order)].copy()

print(
    f"{df['uid'].nunique()} cells | events: {order}"
)


# ==============================================================
# Effect-size helpers
# ==============================================================

def matched_rank_biserial(x_values, y_values):
    """
    Calculate matched-pairs rank-biserial correlation.

    Positive values indicate that y tends to be larger than x.
    """
    differences = (
        np.asarray(y_values, dtype=float)
        - np.asarray(x_values, dtype=float)
    )

    differences = differences[
        np.isfinite(differences)
        & (differences != 0)
    ]

    if differences.size == 0:
        return np.nan

    ranks = (
        pd.Series(np.abs(differences))
        .rank(method="average")
        .to_numpy()
    )

    positive_rank_sum = ranks[
        differences > 0
    ].sum()

    negative_rank_sum = ranks[
        differences < 0
    ].sum()

    total_rank_sum = (
        positive_rank_sum
        + negative_rank_sum
    )

    if total_rank_sum == 0:
        return np.nan

    return (
        positive_rank_sum
        - negative_rank_sum
    ) / total_rank_sum


def kendalls_w(chi_square, n_cells, n_stages):
    """Calculate Kendall's W from the Friedman chi-square statistic."""
    if n_cells <= 0 or n_stages <= 1:
        return np.nan

    return chi_square / (
        n_cells
        * (n_stages - 1)
    )


def safe_paired_wilcoxon(x_values, y_values):
    """Run a paired Wilcoxon test with basic edge-case handling."""
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)

    finite = np.isfinite(x_array) & np.isfinite(y_array)

    x_array = x_array[finite]
    y_array = y_array[finite]

    if len(x_array) < 3:
        return np.nan

    differences = y_array - x_array

    if np.allclose(differences, 0):
        return 1.0

    try:
        return wilcoxon(
            x_array,
            y_array,
            alternative="two-sided",
            zero_method="wilcox",
        ).pvalue

    except ValueError:
        return np.nan


# ==============================================================
# Statistical analysis
# ==============================================================

def run_stats(metric):
    """
    Run Friedman and paired Wilcoxon tests for one metric.

    Pairing is based on the globally unique uid.
    """
    pivot = (
        df.pivot_table(
            index="uid",
            columns="ev",
            values=metric,
            aggfunc="mean",
        )[order]
        .dropna()
    )

    if len(pivot) < 3:
        raise ValueError(
            f"Insufficient complete paired observations for {metric}."
        )

    friedman_chi_square, friedman_p = friedmanchisquare(
        *[
            pivot[event].values
            for event in order
        ]
    )

    kendall_w = kendalls_w(
        friedman_chi_square,
        len(pivot),
        len(order),
    )

    pairs = []
    raw_pvalues = []
    effect_sizes = []

    for first_index in range(len(order)):
        for second_index in range(
            first_index + 1,
            len(order),
        ):
            first_event = order[first_index]
            second_event = order[second_index]

            raw_pvalue = safe_paired_wilcoxon(
                pivot[first_event].values,
                pivot[second_event].values,
            )

            effect_size = matched_rank_biserial(
                pivot[first_event].values,
                pivot[second_event].values,
            )

            pairs.append(
                (
                    first_event,
                    second_event,
                )
            )

            raw_pvalues.append(raw_pvalue)
            effect_sizes.append(effect_size)

    raw_pvalues = np.asarray(
        raw_pvalues,
        dtype=float,
    )

    adjusted_pvalues = np.full(
        raw_pvalues.shape,
        np.nan,
        dtype=float,
    )

    finite = np.isfinite(raw_pvalues)

    if finite.any():
        adjusted_pvalues[finite] = multipletests(
            raw_pvalues[finite],
            method="holm",
        )[1]

    posthoc = list(
        zip(
            pairs,
            raw_pvalues,
            adjusted_pvalues,
            effect_sizes,
        )
    )

    return (
        pivot,
        friedman_chi_square,
        friedman_p,
        kendall_w,
        posthoc,
    )


print(
    "\n=== Friedman + post-hoc paired Wilcoxon with Holm correction ==="
)

print(
    "Omnibus effect size: Kendall's W | "
    "Pairwise effect size: matched rank-biserial r"
)


STATS = {}

for metric, _ in METRICS:
    (
        pivot,
        chi_square,
        p_value,
        kendall_w,
        posthoc,
    ) = run_stats(metric)

    STATS[metric] = (
        pivot,
        chi_square,
        p_value,
        kendall_w,
        posthoc,
    )

    print(
        f"\n{metric}: "
        f"Friedman chi2={chi_square:.2f}, "
        f"p={p_value:.2e}, "
        f"Kendall W={kendall_w:.3f} "
        f"(n={len(pivot)})"
    )

    for (
        first_event,
        second_event,
    ), raw_p, adjusted_p, effect_size in posthoc:

        if not np.isfinite(adjusted_p):
            significance = "n/a"
        elif adjusted_p < 1e-3:
            significance = "***"
        elif adjusted_p < 1e-2:
            significance = "**"
        elif adjusted_p < 0.05:
            significance = "*"
        else:
            significance = "n.s."

        larger_group = (
            second_event
            if effect_size > 0
            else first_event
        )

        print(
            f"    {first_event:<14} vs {second_event:<14}: "
            f"p_raw={raw_p:.2e}, "
            f"p_adj={adjusted_p:.2e} {significance:<4} "
            f"| r={effect_size:+.3f} "
            f"({larger_group} larger, "
            f"|r|={abs(effect_size):.2f})"
        )


# ==============================================================
# Plotting helpers
# ==============================================================

def significance_stars(p_value):
    """Convert an adjusted p-value to significance notation."""
    if not np.isfinite(p_value):
        return "n/a"

    if p_value < 0.001:
        return "***"

    if p_value < 0.01:
        return "**"

    if p_value < 0.05:
        return "*"

    return "n.s."


def style_axis(ax):
    """Apply the common manuscript axis style."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(
        AXIS_LINE_WIDTH
    )

    ax.spines["bottom"].set_linewidth(
        AXIS_LINE_WIDTH
    )

    ax.tick_params(
        axis="x",
        labelsize=X_TICK_LABEL_SIZE,
        width=AXIS_LINE_WIDTH,
        length=3,
        pad=2,
    )

    ax.tick_params(
        axis="y",
        labelsize=Y_TICK_LABEL_SIZE,
        width=AXIS_LINE_WIDTH,
        length=3,
    )


def draw_mean_sd(
    ax,
    position,
    values,
):
    """Draw mean ± SD using an I-bar."""
    mean_value = np.mean(values)
    standard_deviation = np.std(
        values,
        ddof=1,
    )

    # Vertical SD line
    ax.plot(
        [position, position],
        [
            mean_value - standard_deviation,
            mean_value + standard_deviation,
        ],
        color="black",
        linewidth=SD_LINE_WIDTH,
        zorder=4,
        solid_capstyle="butt",
    )

    # Upper SD cap
    ax.plot(
        [
            position - SD_CAP_HALF_WIDTH,
            position + SD_CAP_HALF_WIDTH,
        ],
        [
            mean_value + standard_deviation,
            mean_value + standard_deviation,
        ],
        color="black",
        linewidth=SD_LINE_WIDTH,
        zorder=4,
    )

    # Lower SD cap
    ax.plot(
        [
            position - SD_CAP_HALF_WIDTH,
            position + SD_CAP_HALF_WIDTH,
        ],
        [
            mean_value - standard_deviation,
            mean_value - standard_deviation,
        ],
        color="black",
        linewidth=SD_LINE_WIDTH,
        zorder=4,
    )

    # Mean line
    ax.plot(
        [
            position - MEAN_LINE_HALF_WIDTH,
            position + MEAN_LINE_HALF_WIDTH,
        ],
        [
            mean_value,
            mean_value,
        ],
        color="black",
        linewidth=MEAN_LINE_WIDTH,
        zorder=5,
        solid_capstyle="butt",
    )


def draw_significance_bracket(
    ax,
    x_start,
    x_end,
    y_value,
    bracket_height,
    label,
):
    """Draw one significance bracket."""
    ax.plot(
        [
            x_start,
            x_start,
            x_end,
            x_end,
        ],
        [
            y_value,
            y_value + bracket_height,
            y_value + bracket_height,
            y_value,
        ],
        color="black",
        linewidth=BRACKET_LINE_WIDTH,
        clip_on=False,
        zorder=6,
    )

    ax.text(
        (x_start + x_end) / 2,
        y_value + bracket_height,
        label,
        ha="center",
        va="bottom",
        fontsize=SIGNIFICANCE_FONT_SIZE,
        zorder=7,
    )


def draw_metric_panel(
    ax,
    metric,
    ylabel,
):
    """
    Draw one paired strip-plot panel.

    One jitter value is assigned to each cell and reused across all
    stages, ensuring that each connecting line passes through its points.
    """
    (
        pivot,
        chi_square,
        p_value,
        kendall_w,
        posthoc,
    ) = STATS[metric]

    adjusted_pvalue_lookup = {
        frozenset(pair): adjusted_pvalue
        for (
            pair,
            _,
            adjusted_pvalue,
            _,
        ) in posthoc
    }

    random_generator = np.random.default_rng(
        JITTER_SEED
    )

    n_cells = len(pivot)

    base_positions = np.arange(
        1,
        len(order) + 1,
    )

    cell_jitter = random_generator.uniform(
        -JITTER_WIDTH,
        JITTER_WIDTH,
        size=n_cells,
    )

    all_stage_values = []

    # Draw paired lines first.
    if SHOW_PAIRED_LINES:
        for row_index, (_, row) in enumerate(
            pivot.iterrows()
        ):
            x_line = (
                base_positions
                + cell_jitter[row_index]
            )

            y_line = [
                row[event]
                for event in order
            ]

            ax.plot(
                x_line,
                y_line,
                color=COLOR_PAIRED_LINE,
                linewidth=PAIRED_LINE_WIDTH,
                alpha=PAIRED_LINE_ALPHA,
                zorder=1,
            )

    # Draw scatter points and mean ± SD.
    for event_index, event in enumerate(order):
        values = pivot[event].values
        position = event_index + 1

        all_stage_values.append(values)

        x_scatter = (
            position
            + cell_jitter
        )

        ax.scatter(
            x_scatter,
            values,
            color=COLOR[event],
            alpha=SCATTER_ALPHA,
            s=SCATTER_SIZE,
            edgecolors="none",
            zorder=3,
        )

        draw_mean_sd(
            ax,
            position,
            values,
        )

    y_maximum = max(
        values.max()
        for values in all_stage_values
    )

    y_minimum = min(
        values.min()
        for values in all_stage_values
    )

    y_range = (
        y_maximum
        - y_minimum
    )

    if y_range <= 0:
        y_range = 1.0

    first_bracket_level = (
        y_maximum
        + y_range * 0.08
    )

    bracket_step = (
        y_range
        * 0.13
    )

    bracket_height = (
        y_range
        * 0.02
    )

    # Show adjacent-stage comparisons only.
    for comparison_index in range(
        len(order) - 1
    ):
        first_event = order[
            comparison_index
        ]

        second_event = order[
            comparison_index + 1
        ]

        adjusted_pvalue = adjusted_pvalue_lookup.get(
            frozenset(
                (
                    first_event,
                    second_event,
                )
            ),
            np.nan,
        )

        x_start = comparison_index + 1
        x_end = comparison_index + 2

        y_value = (
            first_bracket_level
            + comparison_index * bracket_step
        )

        draw_significance_bracket(
            ax,
            x_start,
            x_end,
            y_value,
            bracket_height,
            significance_stars(
                adjusted_pvalue
            ),
        )

    # Reserve space for stacked significance brackets.
    y_upper = (
        first_bracket_level
        + (len(order) - 2) * bracket_step
        + bracket_height
        + y_range * 0.12
    )

    ax.set_ylim(
        y_minimum - y_range * 0.08,
        y_upper,
    )

    ax.set_xticks(
        range(
            1,
            len(order) + 1,
        )
    )

    ax.set_xticklabels(
        [
            DISPLAY_LABEL[event]
            for event in order
        ],
        rotation=45,
        ha="center",
        fontsize=X_TICK_LABEL_SIZE,
    )

    ax.set_ylabel(
        ylabel,
        fontsize=AXIS_LABEL_SIZE,
    )

    if SHOW_STATS_IN_TITLE:
        ax.set_title(
            f"{ylabel} "
            f"(W={kendall_w:.2f}, "
            f"p={p_value:.1e})",
            fontsize=TITLE_SIZE,
            pad=4,
        )
    else:
        ax.set_title(
            ylabel,
            fontsize=TITLE_SIZE,
            pad=4,
        )

    style_axis(ax)


# ==============================================================
# Generate Figure 4B-D
# ==============================================================

fig, axes = plt.subplots(
    1,
    3,
    figsize=(
        mm_to_inch(FINAL_FIG_WIDTH_MM),
        mm_to_inch(FINAL_FIG_HEIGHT_MM),
    ),
)

for axis, (metric, ylabel) in zip(
    axes,
    METRICS,
):
    draw_metric_panel(
        axis,
        metric,
        ylabel,
    )

# Manual layout preserves the physical dimensions of the SVG.
fig.subplots_adjust(
    left=0.070,
    right=0.985,
    bottom=0.315,
    top=0.900,
    wspace=0.42,
)


# ==============================================================
# Save optional outputs
# ==============================================================

if SAVE_SVG or SAVE_PDF or SAVE_STATS_CSV:
    os.makedirs(
        SAVE_DIR,
        exist_ok=True,
    )


output_base = os.path.join(
    SAVE_DIR,
    OUTPUT_BASENAME,
)


if SAVE_SVG:
    svg_path = output_base + ".svg"

    fig.savefig(
        svg_path,
        transparent=True,
    )

    print(f"\nSaved SVG: {svg_path}")


if SAVE_PDF:
    pdf_path = output_base + ".pdf"

    fig.savefig(
        pdf_path,
        transparent=True,
    )

    print(f"Saved PDF: {pdf_path}")


# ==============================================================
# Descriptive statistics and optional CSV
# ==============================================================

statistics_records = []

print(
    "\n=== Per-stage descriptive statistics for Figure 4 ==="
)

for metric, _ in METRICS:
    (
        pivot,
        chi_square,
        p_value,
        kendall_w,
        posthoc,
    ) = STATS[metric]

    print(
        f"\n{metric} "
        f"(n={len(pivot)} paired cells, "
        f"pooled across {n_experiments} experiments | "
        f"Kendall W={kendall_w:.3f}):"
    )

    for event in order:
        values = pivot[event].values

        mean_value = values.mean()
        sd_value = values.std(ddof=1)
        median_value = np.median(values)

        print(
            f"    {event:<14}: "
            f"mean {mean_value:.3f} ± {sd_value:.3f} "
            f"| median {median_value:.3f}"
        )

        statistics_records.append({
            "metric": metric,
            "result_type": "descriptive",
            "comparison_or_stage": event,
            "n": len(values),
            "mean": mean_value,
            "sd": sd_value,
            "median": median_value,
            "friedman_chi_square": chi_square,
            "friedman_p": p_value,
            "kendall_w": kendall_w,
            "raw_p": np.nan,
            "adjusted_p": np.nan,
            "rank_biserial_r": np.nan,
        })

    for (
        first_event,
        second_event,
    ), raw_p, adjusted_p, effect_size in posthoc:

        statistics_records.append({
            "metric": metric,
            "result_type": "pairwise",
            "comparison_or_stage": (
                f"{first_event} vs {second_event}"
            ),
            "n": len(pivot),
            "mean": np.nan,
            "sd": np.nan,
            "median": np.nan,
            "friedman_chi_square": chi_square,
            "friedman_p": p_value,
            "kendall_w": kendall_w,
            "raw_p": raw_p,
            "adjusted_p": adjusted_p,
            "rank_biserial_r": effect_size,
        })

        print(
            f"      {first_event} vs {second_event}: "
            f"p_adj={adjusted_p:.3e}, "
            f"r={effect_size:+.3f}"
        )


if SAVE_STATS_CSV:
    statistics_path = (
        output_base
        + "_statistics.csv"
    )

    pd.DataFrame(
        statistics_records
    ).to_csv(
        statistics_path,
        index=False,
    )

    print(
        f"\nSaved statistics CSV: "
        f"{statistics_path}"
    )


plt.show()
