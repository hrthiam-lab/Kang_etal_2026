# NucMetrics

An ImageJ/Fiji macro toolset for computing DNA staining-based image metrics for live-cell tracking of chromatin organization.

NucMetrics computes three metrics from nuclear DNA staining images:
- **CV** (Coefficient of Variation) — computed from raw pixel intensities
- **1-Gini** (complement of Gini coefficient) — computed from min-max normalized intensities
- **DSI** (Diffuse Signal Index) — computed from min-max normalized intensities

## Installation

1. Download `NucMetrics_Toolset.ijm`
2. Copy it to your Fiji (or ImageJ) installation: `Fiji.app/macros/toolsets/`
3. Restart Fiji/ImageJ
4. Click `>>` on the toolbar and select **NucMetrics_Toolset**
5. The NucMetrics icon appears on the toolbar — click to run

**Requirements:** Fiji or ImageJ (1.53b or later). No additional plugins or dependencies.

## Modes

### Mode 1: Current Selection

Draw a freehand, polygon, or oval ROI around a nucleus, then run NucMetrics. If no ROI is drawn, NucMetrics will activate the Freehand tool and prompt you to draw one.

### Mode 2: ROI Manager

Pre-load multiple ROIs into the ROI Manager (must be open before running), then run NucMetrics. All ROIs are batch-processed and results are output to the Results table.

> **Note:** The ROI Manager window must be open with ROIs added before selecting this mode.

### Mode 3: Binary Mask

Use an external binary mask image to define nuclear regions. The mask image must already be open in Fiji/ImageJ alongside the original DNA-stained image.

- **Single slice:** Select a mask image from the dropdown and compute metrics for the current slice.
- **Whole stack:** If both the original image and the mask are stacks with matching dimensions, NucMetrics can process all corresponding planes. Each mask plane defines the nuclear region for the matching plane in the original image. This mode measures one nucleus per plane and is intended for single-nucleus stacks.

> **Note:** The mask image must be opened in Fiji/ImageJ before running NucMetrics. "Binary mask" refers to selecting an already-open mask window, not importing a file.

### Mode 4: Auto-Generate Binary Mask

Automatically generate a binary mask using intensity thresholding (Li, Otsu, or Triangle), followed by morphological cleanup (fill holes, opening) and particle analysis.

**Workflow:**
1. NucMetrics generates a binary mask and displays it in a new window
2. Review the mask visually
3. Click OK to compute metrics using the generated mask, or Cancel to keep the mask only for manual editing

**Stack support:**
- **Single slice:** Generates a mask for the currently displayed slice
- **Entire series:** Generates a mask stack across all time points (T) or Z-slices. For hyperstacks with both T and Z dimensions, you can choose which axis to iterate over while keeping the other fixed at the current position
- Whole-stack mode is designed for **single-nucleus** stacks (one nucleus per field of view). It does not perform multi-object tracking across time points

> **Tip:** For multi-nuclei time-lapse data, use Mode 4 to generate a binary mask stack, manually review/correct it, then use Mode 3 (Binary Mask, whole stack) to compute metrics.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| DSI threshold (tau) | 0.3 | Normalized intensity threshold for DSI computation |
| Min nucleus area | 30 px | Minimum object size for segmentation modes |
| Auto-threshold method | Li | Thresholding algorithm (Li, Otsu, or Triangle) |

## How metrics are computed

All metrics are computed from pixel intensities within the nuclear ROI.

**CV** = sigma / mu, where sigma and mu are the standard deviation and mean of raw pixel intensities.

**1-Gini**: Raw intensities are min-max normalized to [0, 1]. The Gini coefficient is computed as:

```
Gini = [2 * sum(i * x_sorted_i)] / [N * sum(x)] - (N + 1) / N
```

The reported value is 1 - Gini, bounded between 0 (maximally unequal) and 1 (perfectly uniform).

**DSI**: Raw intensities are min-max normalized to [0, 1]. DSI is the fraction of normalized pixels exceeding the threshold tau:

```
DSI = count(x_i > tau) / N
```

DSI ranges from 0 to 1, with higher values indicating more spatially uniform DNA signal distribution.

## Output

Results are displayed in the ImageJ Results table with the following columns:

| Column | Description |
|--------|-------------|
| Label | Nucleus identifier (e.g., Auto_1, T001_N1) |
| N_pixels | Number of pixels in the nuclear ROI |
| Mean_raw | Mean raw pixel intensity |
| StdDev_raw | Standard deviation of raw pixel intensity |
| CV | Coefficient of variation |
| 1-Gini | Complement of the Gini coefficient |
| DSI | Diffuse Signal Index |

For whole-stack processing, labels follow the format `T###_N1` (time point) or `Z###_N1` (Z-slice) or `S###_N1` (plain stack slice).

Results can be exported via `File > Save As...` from the Results table.

## Citation

If you use NucMetrics, please cite:

> Kang M, Cabral AT, Sawant M, Thiam HR (2026). Benchmarking three simple DNA staining-based image metrics for live-cell tracking of chromatin organization. *bioRxiv*. [https://doi.org/10.64898/2026.03.30.715467](https://doi.org/10.64898/2026.03.30.715467)

## License

MIT
