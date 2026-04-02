# Kang_2026_DSI_paper

Code for figure generation, CV/1-Gini/DSI analysis, and the **NucMetrics** ImageJ/Fiji toolset accompanying:

> Kang M, Cabral AT, Sawant M, Thiam HR (2026). Benchmarking three simple DNA staining-based image metrics for live-cell tracking of chromatin organization. *bioRxiv*. [https://doi.org/10.64898/2026.03.30.715467](https://doi.org/10.64898/2026.03.30.715467)

## Overview

This repository contains:

- **NucMetrics** — an ImageJ/Fiji macro toolset for computing CV, 1-Gini, and DSI on DNA-stained nuclear images
- **Figure scripts** — Python scripts for reproducing all figures in the manuscript
- **Example data** — sample images for testing

## NucMetrics

NucMetrics is a single-file ImageJ/Fiji macro toolset that computes three DNA staining-based image metrics directly within Fiji, with no programming required.

**Metrics:**
- **CV** (Coefficient of Variation) — computed from raw pixel intensities
- **1-Gini** (complement of Gini coefficient) — computed from min-max normalized intensities
- **DSI** (Diffuse Signal Index) — computed from min-max normalized intensities

**Modes:**

| Mode | Description |
|------|-------------|
| Current Selection | Draw a freehand or polygon ROI on a nucleus |
| ROI Manager | Batch-process multiple pre-defined ROIs |
| Binary Mask | Import an external segmentation mask |
| Auto-Segmentation | Automatic thresholding (Li, Otsu, or Triangle) |

### Installation

1. Download [`NucMetrics/NucMetrics_Toolset.ijm`](NucMetrics/NucMetrics_Toolset.ijm)
2. Copy to `Fiji.app/macros/toolsets/`
3. Restart Fiji
4. Click `>>` on the toolbar and select **NucMetrics_Toolset**

**Requirements:** Fiji (ImageJ 1.53b or later). No additional plugins or dependencies.

### Usage

1. Open a DNA-stained image in Fiji
2. Click the NucMetrics icon on the toolbar
3. Select a mode and set the DSI threshold (default: τ = 0.3)
4. Results appear in the ImageJ Results table

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| DSI threshold (τ) | 0.3 | Normalized intensity threshold for DSI |
| Min nucleus area | 30 px | Minimum object size for auto-segmentation |
| Auto-threshold method | Li | Thresholding algorithm (Li, Otsu, or Triangle) |

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
│   └── README.md                 # Detailed NucMetrics documentation
├── Figures/
│   ├── fig1_metric_explanation.py
│   ├── fig2_live_cell_tracking.py
│   ├── fig3_tn5_correlation.py
│   └── suppfig1_threshold_sweep.py
├── example_data/                 # Sample images for testing
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
