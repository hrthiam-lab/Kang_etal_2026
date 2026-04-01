# NucMetrics

An ImageJ/Fiji macro toolset for computing DNA staining-based image metrics for live-cell tracking of chromatin organization.

NucMetrics computes three metrics from nuclear DNA staining images:
- **CV** (Coefficient of Variation) — computed from raw pixel intensities
- **1-Gini** (complement of Gini coefficient) — computed from min-max normalized intensities
- **DSI** (Diffuse Signal Index) — computed from min-max normalized intensities

## Installation

1. Download `NucMetrics_Toolset.ijm`
2. Copy it to your Fiji installation: `Fiji.app/macros/toolsets/`
3. Restart Fiji
4. Click `>>` on the toolbar and select **NucMetrics_Toolset**
5. The NucMetrics icon appears on the toolbar — click to run

**Requirements:** Fiji (ImageJ 1.53b or later). No additional plugins or dependencies.

## Modes

| Mode | Description |
|------|-------------|
| **Current Selection** | Draw a freehand or polygon ROI on a nucleus, then run |
| **ROI Manager** | Pre-load multiple ROIs into ROI Manager for batch processing |
| **Binary Mask** | Open a binary mask image alongside the original; nuclear regions are extracted automatically |
| **Auto-Segmentation** | Automatic thresholding (Li, Otsu, or Triangle) with morphological cleanup and particle analysis |

## Usage

1. Open a DNA-stained image in Fiji
2. Click the NucMetrics icon on the toolbar
3. Select a mode and set the DSI threshold (default: τ = 0.3)
4. For Auto-Segmentation: a second dialog allows selection of threshold method and minimum nucleus area
5. Results (CV, 1-Gini, DSI, mean intensity, standard deviation, pixel count) appear in the Results table
6. Export results via `File > Save As...` from the Results table

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| DSI threshold (τ) | 0.3 | Normalized intensity threshold for DSI computation |
| Min nucleus area | 30 px | Minimum object size for auto-segmentation and binary mask modes |
| Auto-threshold method | Li | Thresholding algorithm (Li, Otsu, or Triangle) |

## How metrics are computed

All metrics are computed from pixel intensities within the nuclear ROI.

**CV** = σ / μ, where σ and μ are the standard deviation and mean of raw pixel intensities.

**1-Gini**: Raw intensities are min-max normalized to [0, 1]. The Gini coefficient is computed as:

```
Gini = [2 * Σ(i * x_sorted_i)] / [N * Σ(x)] - (N + 1) / N
```

The reported value is 1 − Gini, bounded between 0 (maximally unequal) and 1 (perfectly uniform).

**DSI**: Raw intensities are min-max normalized to [0, 1]. DSI is the fraction of normalized pixels exceeding the threshold τ:

```
DSI = count(x_i > τ) / N
```

DSI ranges from 0 to 1, with higher values indicating more spatially uniform DNA signal distribution.

## Citation

If you use NucMetrics, please cite:

> Kang M, Cabral AT, Sawant M, Thiam HR (2026). Benchmarking three simple DNA staining-based image metrics for live-cell tracking of chromatin organization. *bioRxiv*. [https://doi.org/10.64898/2026.03.30.715467](https://doi.org/10.64898/2026.03.30.715467)

## License

MIT
