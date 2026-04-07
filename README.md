# Kang et al., 2026

Code for figure generation, CV/1-Gini/DSI analysis, and the **NucMetrics** ImageJ/Fiji toolset accompanying:

> Kang M, Cabral AT, Sawant M, Thiam HR (2026). Benchmarking three simple DNA staining-based image metrics for live-cell tracking of chromatin organization. *bioRxiv*. [https://doi.org/10.64898/2026.03.30.715467](https://doi.org/10.64898/2026.03.30.715467)

## Overview

This repository contains:

- **NucMetrics** — an ImageJ/Fiji macro toolset for computing CV, 1-Gini, and DSI on DNA-stained nuclear images
- **Figure scripts** — Python scripts for reproducing all figures in the manuscript
- **Example data** — sample images and masks for testing NucMetrics and validating outputs against the Python analysis code

## NucMetrics

NucMetrics is a single-file ImageJ/Fiji macro toolset that computes three DNA
staining-based image metrics that quantify intranuclear signal distribution as
a readout of chromatin organization, directly within Fiji with no programming required.

![NucMetrics toolbar icon](NucMetrics/NucMetrics_toolbar.png)

**Metrics:**
- **CV** (Coefficient of Variation) — from raw pixel intensities
- **1-Gini** (complement of Gini coefficient) — from min-max normalized intensities
- **DSI** (Diffuse Signal Index) — from min-max normalized intensities

**Modes:**
> ⚠️ Whole-stack mode supports single-nucleus stacks only. See [NucMetrics/README.md](NucMetrics/README.md) for details on multi-nuclei data.

| Mode | Description |
|------|-------------|
| Current Selection | Draw an ROI on a nucleus manually |
| ROI Manager | Batch-process multiple pre-defined ROIs |
| Binary Mask | Use an external binary mask image (single slice or matched stack) |
| Auto-Generate Binary Mask | Automatic thresholding, mask review, then compute (single slice or entire series) |

### Quick start

There are two ways to install NucMetrics. See [NucMetrics/README.md](NucMetrics/README.md) for full documentation and a tutorial.

**Method A — Install to toolsets folder (recommended for repeated use):**

1. Download [`NucMetrics/NucMetrics_Toolset.ijm`](NucMetrics/NucMetrics_Toolset.ijm)
2. Save it in the Fiji application folder at: **`Fiji.app/macros/toolsets/NucMetrics_Toolset.ijm`**
 
   > **💡 How to find the `toolsets` folder:**
   > - **macOS:** In Finder, locate `Fiji.app` (check your `Downloads`, `Applications`, or wherever you saved it). Right-click → **Show Package Contents** → `macros` → `toolsets`. If Fiji is in your Dock, right-click the Dock icon → `Options` → `Show in Finder`.
   > - **Windows:** Press the Start key and type `fiji`. Right-click the result → **Open file location** to open the `Fiji.app` folder → `macros` → `toolsets`
   > - **Linux:** Open the `Fiji.app` directory → `macros` → `toolsets`
   >
   > ⚠️ The file must be inside the **`toolsets`** subfolder, NOT directly in `macros`. This is also **not** the same as `Plugins > Macros` in the Fiji menu, which points to a different folder.

3. Click `>>` on the toolbar and select **NucMetrics_Toolset** — no restart required
4. Open a DNA-stained image, then click the NucMetrics icon to run

**Method B — Temporary install (current session only):**

1. In Fiji, go to `Plugins > Macros > Install...` and select the downloaded `.ijm` file
2. The NucMetrics icon appears on the toolbar for the current session

**Requirements:** Fiji or ImageJ (1.53b or later). No additional plugins or dependencies.

## Figure scripts

Python scripts for generating all manuscript figures are in the [`Figures/`](Figures/) folder.

| Script | Figure |
|--------|--------|
| `fig1_metric_explanation.py` | Fig. 1 — Metric definitions and illustration |
| `fig2_live_cell_tracking.py` | Fig. 2 — Live-cell trajectory analysis |
| `fig3_tn5_correlation.py` | Fig. 3 — Correlation with Tn5/ATAC-see |
| `suppfig1_threshold_sweep.py` | Supp. Fig. 1 — DSI threshold selection |

## Repository structure

```
Kang_2026_DSI_paper/
├── NucMetrics/
│   ├── NucMetrics_Toolset.ijm    # Fiji macro toolset
│   ├── NucMetrics_toolbar.png    # Toolbar icon screenshot
│   └── README.md                 # Detailed documentation and tutorial
├── Figures/
│   ├── fig1_metric_explanation.py
│   ├── fig2_live_cell_tracking.py
│   ├── fig3_tn5_correlation.py
│   └── suppfig1_threshold_sweep.py
├── example_data/                 # Sample images and masks for testing and validation
├── LICENSE
└── README.md
```

## Citation

If you use NucMetrics or any code from this repository, please cite:

```bibtex
@article{Kang2026,
  author  = {Kang, Minwoo and Cabral, Aidan Tomas and Sawant, Manasi and Thiam, Hawa Racine},
  title   = {Benchmarking three simple DNA staining-based image metrics for live-cell tracking of chromatin organization},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.03.30.715467}
}
```

## License

[MIT](LICENSE)
