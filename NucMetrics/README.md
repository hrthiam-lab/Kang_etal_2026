# NucMetrics

An ImageJ/Fiji macro toolset for computing DNA staining-based image metrics for live-cell tracking of chromatin organization.

NucMetrics computes three metrics from nuclear DNA staining images:
- **CV** (Coefficient of Variation) — computed from raw pixel intensities
- **1-Gini** (complement of Gini coefficient) — computed from min-max normalized intensities
- **DSI** (Diffuse Signal Index) — computed from min-max normalized intensities

## Installation

There are two ways to install NucMetrics. **Method A** is recommended for repeated use.

### Method A: Install to the toolsets folder (recommended)

This makes NucMetrics permanently available in the `>>` toolbar menu.

1. Download `NucMetrics_Toolset.ijm`
2. Save it in the Fiji application folder at:
 
   **`Fiji.app/macros/toolsets/NucMetrics_Toolset.ijm`**
 
   > **💡 How to find the `toolsets` folder:**
   > - **macOS:** In Finder, locate `Fiji.app` (check your `Downloads`, `Applications`, or wherever you saved it). Right-click → **Show Package Contents** → `macros` → `toolsets`. If Fiji is in your Dock, right-click the Dock icon → `Options` → `Show in Finder`.
   > - **Windows:** Press the Start key and type `fiji`. Right-click the result → **Open file location** to open the `Fiji.app` folder → `macros` → `toolsets`
   > - **Linux:** Open the `Fiji.app` directory → `macros` → `toolsets`
   >
   > ⚠️ The file must be inside the **`toolsets`** subfolder, NOT directly in `macros`. This is also **not** the same as `Plugins > Macros` in the Fiji menu, which points to a different folder (`plugins/Macros`).
 
3. Open Fiji (no restart needed if already open). Click the **`>>`** button on the right end of the toolbar and select **NucMetrics_Toolset**.

4. The NucMetrics icon appears on the toolbar:

   ![NucMetrics toolbar icon](NucMetrics_toolbar.png)

> **Note:** Each time you open Fiji, you need to click `>>` and select NucMetrics_Toolset to load it into the toolbar. The same procedure can be used for ImageJ.

### Method B: Temporary install via menu (quick start)

This loads NucMetrics for the current session only. You will need to repeat this each time you open Fiji.

1. Download `NucMetrics_Toolset.ijm` to any location on your computer
2. In Fiji, go to `Plugins > Macros > Install...`
3. Navigate to and select `NucMetrics_Toolset.ijm`
4. The NucMetrics icon appears on the toolbar for the current session

### First use

> ⚠️ **Important:** Always open a DNA-stained image **before** clicking the NucMetrics icon. If no image is open, NucMetrics will show an error.

1. Open a DNA-stained image in Fiji (`File > Open...`)
2. Click the NucMetrics icon on the toolbar (or press `1` on the numpad)
3. Select a mode, set the DSI threshold (tau), and click OK

**Requirements:** Fiji or ImageJ (1.53b or later). No additional plugins or dependencies.

## Tutorial: Quick test with example data

Try NucMetrics using the example images provided in this repository.

### Test 1: Single nucleus with Mode 1 (Current Selection)

1. Open `example_data/compact.tif` in Fiji
2. Click the NucMetrics icon
3. Select **Current Selection** mode, keep tau = 0.3, click OK
4. NucMetrics will prompt you to draw an ROI — draw a freehand selection around the nucleus, then click OK
5. Check the Results table. Expected values (approximate):
   - CV ≈ 0.58, 1-Gini ≈ 0.62, DSI ≈ 0.37

### Test 2: Single nucleus with Mode 3 (Binary Mask)

1. Open `example_data/compact.tif` in Fiji
2. Also open `example_data/compact_mask.tif`
3. Click on the 'compact.tif' image window to make it active
4. Click the NucMetrics icon
5. Select **Binary Mask** mode → select the mask image from the dropdown (compact_mask.tif)→ click OK
6. Check the Results table for CV, 1-Gini, and DSI values

### Test 3: Multiple nuclei with Mode 4 (Auto-Generate Binary Mask)

1. Open `example_data/multiple_nuc.tif` in Fiji
2. Click the NucMetrics icon
3. Select **Auto-Generate Binary Mask** mode → choose a threshold method (Li is the default and works well for this example; try Li, Otsu, or Triangle and compare) → click OK
4. Review the generated mask in the new window. The review dialog does not block the image windows, so you can zoom, pan, and compare the mask against the original while it is open
5. In the review dialog, choose one of:
   - **Compute NucMetrics with this mask** — measure using this mask
   - **Go back and adjust the mask settings** — return to the settings dialog, change the threshold method / min area / edge exclusion, and regenerate the mask (no need to cancel and start over)
   - **Stop and keep the generated mask only** — keep the mask window for manual editing

> **Tip:** The built-in auto-thresholding may not work well for all datasets. If you have a segmentation pipeline that produces better masks for your data, save the mask as a binary image and use **Mode 3 (Binary Mask)** instead.

> **Validation:** The example images (`compact.tif`, `decompact.tif`) are the same ones used in Fig. 1 of the manuscript. NucMetrics outputs on these images have been validated against the Python analysis code to produce identical CV, 1-Gini, and DSI values.

### Test 4: Time-lapse trajectory with Mode 4 (whole-stack processing)

This test demonstrates frame-by-frame batch processing of an entire time-lapse movie for trajectory quantification.

1. Open `example_data/Single_nuc_timelapse.tif` in Fiji (a 49-frame time-lapse of a single dHL-60 nucleus undergoing NETosis)
2. Click the NucMetrics icon
3. Select **Auto-Generate Binary Mask** mode → click OK
4. In the settings dialog:
   - **Threshold method:** Li (default)
   - **Min nucleus area:** 30 (default)
   - Under **Stack options**, check **"Generate mask for an entire series"**
   - Click OK
5. NucMetrics generates a binary mask stack for all 49 frames and displays it in a new window
6. The "Review Auto Mask" dialog opens without blocking the image windows. Scroll through all 49 mask planes with the stack slider, the `<` / `>` keys, or the mouse wheel to verify segmentation quality across time points
7. In the review dialog, select **Compute NucMetrics with this mask** and click OK. If the segmentation looks wrong, select **Go back and adjust the mask settings** instead to change the threshold method or minimum area and regenerate the stack
8. The Results table now contains 49 rows, labeled `T001_N1` through `T049_N1`, with CV, 1-Gini, and DSI for each frame
9. Export the Results table (`File > Save As...`) for trajectory plotting in Python, R, or Excel

> **Note:** This is the recommended workflow for trajectory-level analysis. The same approach works with Mode 3 (Binary Mask) if you have pre-made mask stacks from an external segmentation pipeline.

## Trajectory quantification workflow

NucMetrics supports frame-by-frame batch processing of time-lapse stacks for trajectory-level analysis. This is useful for tracking how CV, 1-Gini, and DSI change over time in a single nucleus (e.g., during NETosis, mitotic exit, or drug treatment).

### Step-by-step: Time-lapse trajectory analysis

1. **Prepare a single-nucleus stack.** Crop your time-lapse movie so that each stack contains one nucleus across all time points.

2. **Open the stack in Fiji** (`File > Open...`)

3. **Choose a processing mode:**

   **Option A — Auto-generate masks (Mode 4):**
   - Click the NucMetrics icon → select **Auto-Generate Binary Mask**
   - Select **Entire series** when prompted
   - Scroll through the generated mask stack in the non-blocking review dialog → choose **Compute NucMetrics with this mask**, or **Go back and adjust the mask settings** to retry with different parameters

   **Option B — Use pre-made masks (Mode 3):**
   - Open your binary mask stack alongside the image stack (must have matching dimensions)
   - Click the NucMetrics icon → select **Binary Mask**
   - Select the mask from the dropdown → check **Whole stack** → click OK

4. **Read the Results table.** Each time point appears as a row labeled `T001_N1`, `T002_N1`, `T003_N1`, etc.

5. **Export for downstream analysis.** Copy or save (`File > Save As...`) the Results table for trajectory plotting in Python, R, MATLAB, or Excel.

> **Note:** Whole-stack mode processes one nucleus per frame. For multi-nuclei fields of view, crop individual nuclei into separate stacks first, or generate per-nucleus mask stacks using an external segmentation/tracking pipeline and use Mode 3.

## Modes

> ⚠️ **Stack processing note:** For single time-point images, multiple nuclei per field of view are fully supported. For time-series (stack) data, the field of view should contain a **single nucleus** — whole-stack mode does not perform multi-object tracking across time points. For multi-nuclei time-lapse data, either (1) crop individual nuclei into separate single-nucleus stacks, or (2) generate per-nucleus binary mask stacks using your own segmentation/tracking pipeline, then use Mode 3 (Binary Mask, whole stack) to compute metrics.

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
1. Set the mask parameters (threshold method, minimum nucleus area, edge exclusion, stack options)
2. NucMetrics generates a binary mask and displays it in a new window
3. Review the mask. The review dialog is non-blocking: the image windows stay interactive, so a generated mask **stack** can be scrolled plane by plane and compared with the original before you commit to it
4. Choose what happens next:

| Review option | Effect |
|---------------|--------|
| Compute NucMetrics with this mask | Measures CV, 1-Gini and DSI using the generated mask |
| Go back and adjust the mask settings | Reopens the settings dialog, regenerates the mask, and returns to this review step — repeat as often as needed |
| Stop and keep the generated mask only | Leaves the mask window open for manual editing (then use Mode 3) |

If the generated mask contains no objects at all, the review dialog says so up front and preselects the "adjust settings" option.

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
