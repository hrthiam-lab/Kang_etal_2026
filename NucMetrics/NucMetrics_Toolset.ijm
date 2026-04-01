// =========================================================================
// NucMetrics_Toolset.ijm
//
// NucMetrics: An ImageJ/Fiji macro toolset for computing DNA staining-
// based image metrics for live-cell tracking of chromatin organization.
//
// Metrics computed:
//   1. CV  (Coefficient of Variation) = sigma / mu  (raw intensities)
//   2. 1-Gini (complement of Gini coefficient, min-max normalized)
//   3. DSI (DNA Staining Index, min-max normalized)
//
// Modes:
//   1. Current Selection  - compute on the active ROI
//   2. ROI Manager        - batch compute for all ROIs in ROI Manager
//   3. Binary Mask        - use an open binary mask image to define nuclei
//   4. Auto-Segmentation  - threshold + Analyze Particles
//
// Installation:
//   1. Copy this file to: Fiji.app/macros/toolsets/
//   2. Restart Fiji
//   3. Click ">>" on the toolbar, select "NucMetrics_Toolset"
//   4. Click the NM icon to run
//   5. Double-click the NM icon to change default settings
//
// If you use this tool, please cite:
//   Kang M, Cabral AT, Sawant M, Thiam HR (2026).
//   "Benchmarking three simple DNA staining-based image metrics for
//   live-cell tracking of chromatin organization."
//   bioRxiv. https://doi.org/10.64898/2026.03.30.715467
//
// Created: 2026-04-01
// Authors: Minwoo Kang Ph.D.
//          HR Thiam Lab, Stanford University
// License: MIT
// =========================================================================

// ---- Global parameters (editable via double-click on icon) ----
var DSI_THRESHOLD = 0.3;
var MIN_AREA = 30;
var AUTO_THRESH_METHOD = "Li";

// =========================================================================
// TOOLBAR ICON: "Nuc" (top) / "Met" (bottom) in blue
// T = text command: Txy ss "char"
//   C05f = dark blue
//   Top row  "Nuc": N at (0,6), u at (5,6), c at (9,6)  size=07
//   Bot row  "Met": M at (0,e), e at (5,e), t at (9,e)  size=07
// =========================================================================
macro "NucMetrics Action Tool - C05fD03D04D05D06D07D08D09D0aD0bD0cD12D13D1cD1dD21D22D2dD2eD30D31D3eD3fD40D4fD50D5fD60D6fD70D7fD80D8fD90D9fDa0DafDb0Db1DbeDbeDbfDc1Dc2DcdDceDd2Dd3DdcDddDe3De4De5De6De7De8De9DeaDebDecCadfD14D15D16D17D18D19D1aD1bD23D24D25D26D27D28D29D2aD2bD2cD32D33D34D35D36D37D38D39D3aD3bD3cD3dD41D42D43D44D45D46D47D48D49D4aD4bD4cD4dD4eD51D52D53D54D55D56D57D58D59D5aD5bD5cD5dD5eD61D62D63D64D65D66D67D68D69D6aD6bD6cD6dD6eD71D72D73D74D75D76D77D78D79D7aD7bD7cD7dD7eD81D82D83D84D85D86D87D88D89D8aD8bD8cD8dD8eD91D92D93D94D95D96D97D98D99D9aD9bD9cD9dD9eDa1Da2Da3Da4Da5Da6Da7Da8Da9DaaDabDacDadDaeDb2Db3Db4Db5Db6Db7Db8Db9DbaDbbDbcDbdDc3Dc4Dc5Dc6Dc7Dc8Dc9DcaDcbDccDd4Dd5Dd6Dd7Dd8Dd9DdaDdbC005T3a08nT8a08uTda08c" {
    runNucMetrics();
}

macro "NucMetrics Action Tool Options" {
    Dialog.create("NucMetrics - Settings");
    Dialog.addMessage("Default parameters for NucMetrics.\n"
        + "These can also be changed at run time.");
    Dialog.addNumber("DSI threshold (tau):", DSI_THRESHOLD);
    Dialog.addNumber("Min nucleus area (px):", MIN_AREA);
    Dialog.addChoice("Auto-threshold method:",
        newArray("Li", "Otsu", "Triangle"), AUTO_THRESH_METHOD);
    Dialog.show();
    DSI_THRESHOLD = Dialog.getNumber();
    MIN_AREA = Dialog.getNumber();
    AUTO_THRESH_METHOD = Dialog.getChoice();
}

// Keyboard shortcut: numpad 1
macro "NucMetrics [n1]" {
    runNucMetrics();
}

// =========================================================================
// MAIN: Two-step dialog - mode first, then mode-specific options
// =========================================================================
function runNucMetrics() {

    if (nImages == 0) {
        exit("No image open.\nPlease open a DNA-stained image first.");
    }

    originalID = getImageID();

    // Save current slice position for stacks
    origSlice = 1;
    origChannel = 1;
    origFrame = 1;
    if (nSlices > 1) {
        Stack.getPosition(origChannel, origSlice, origFrame);
    }

    // --- STEP 1: Mode + DSI threshold ---
    modes = newArray(
        "Current Selection",
        "ROI Manager",
        "Binary Mask",
        "Auto-Segmentation"
    );

    Dialog.create("NucMetrics");
    Dialog.addMessage("NucMetrics\n"
        + "Compute CV, 1-Gini, and DSI\n"
        + "on DNA-stained nuclear images.\n ");
    Dialog.addChoice("Mode:", modes, modes[0]);
    Dialog.addNumber("DSI threshold (tau):", DSI_THRESHOLD);
    Dialog.addHelp("https://doi.org/10.64898/2026.03.30.715467");
    Dialog.show();

    mode = Dialog.getChoice();
    DSI_THRESHOLD = Dialog.getNumber();

    // --- STEP 2: Mode-specific dialog ---
    if (mode == "Auto-Segmentation") {
        Dialog.create("NucMetrics - Auto-Segmentation");
        Dialog.addChoice("Threshold method:",
            newArray("Li", "Otsu", "Triangle"), AUTO_THRESH_METHOD);
        Dialog.addNumber("Min nucleus area (px):", MIN_AREA);
        Dialog.addCheckbox("Exclude nuclei touching edges", true);
        Dialog.show();
        AUTO_THRESH_METHOD = Dialog.getChoice();
        MIN_AREA = Dialog.getNumber();
        excludeEdges = Dialog.getCheckbox();
        modeAutoSegmentation(originalID, excludeEdges,
            origChannel, origSlice, origFrame);

    } else if (mode == "Binary Mask") {
        maskList = newArray(0);
        maskIDs = newArray(0);
        for (i = 1; i <= nImages; i++) {
            selectImage(i);
            t = getTitle();
            id = getImageID();
            if (id != originalID) {
                maskList = Array.concat(maskList, t);
                maskIDs = Array.concat(maskIDs, id);
            }
        }
        selectImage(originalID);

        if (maskList.length == 0) {
            exit("No other image is open.\n"
                + "Open a binary mask image, then try again.");
        }

        Dialog.create("NucMetrics - Binary Mask");
        Dialog.addChoice("Mask image:", maskList, maskList[0]);
        Dialog.addNumber("Min nucleus area (px):", MIN_AREA);
        Dialog.show();
        maskChoice = Dialog.getChoice();
        MIN_AREA = Dialog.getNumber();

        maskID = 0;
        for (i = 0; i < maskList.length; i++) {
            if (maskList[i] == maskChoice) {
                maskID = maskIDs[i];
            }
        }
        modeBinaryMask(originalID, maskID);

    } else if (mode == "ROI Manager") {
        modeROIManager(originalID);

    } else {
        modeCurrentSelection(originalID);
    }
}

// =========================================================================
// MODE 1: Current Selection
// =========================================================================
function modeCurrentSelection(imgID) {
    selectImage(imgID);
    if (selectionType() == -1) {
        exit("No selection found.\nDraw an ROI on a nucleus first.");
    }
    pixels = getPixelsFromSelection(imgID);
    if (pixels.length == 0) {
        exit("No pixels found in selection.");
    }
    results = computeMetrics(pixels);
    printResults("Current Selection", results);
    run("Clear Results");
    addToResultsTable("Selection", results);
}

// =========================================================================
// MODE 2: ROI Manager
// =========================================================================
function modeROIManager(imgID) {
    n = roiManager("count");
    if (n == 0) {
        exit("ROI Manager is empty.\nAdd nuclear ROIs first.");
    }

    run("Clear Results");
    for (i = 0; i < n; i++) {
        selectImage(imgID);
        roiManager("select", i);
        roiName = Roi.getName();
        if (roiName == "" || roiName == "null") roiName = "ROI_" + (i + 1);

        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            addToResultsTable(roiName, results);
        }
    }
    roiManager("deselect");
    showCompletionMessage(n, "ROI Manager");
}

// =========================================================================
// MODE 3: Binary Mask
// =========================================================================
function modeBinaryMask(imgID, maskID) {
    selectImage(imgID);
    w1 = getWidth(); h1 = getHeight();
    selectImage(maskID);
    w2 = getWidth(); h2 = getHeight();

    if (w1 != w2 || h1 != h2) {
        exit("Dimension mismatch.\n"
            + "Image: " + w1 + "x" + h1 + "\n"
            + "Mask:  " + w2 + "x" + h2);
    }

    // Extract ROIs from mask (handle stacks)
    selectImage(maskID);
    if (nSlices > 1) {
        Stack.getPosition(ch, sl, fr);
        run("Duplicate...", "title=_NM_temp_mask_ duplicate channels=" + ch
            + " slices=" + sl + " frames=" + fr);
    } else {
        run("Duplicate...", "title=_NM_temp_mask_");
    }
    tempMaskID = getImageID();
    while (nSlices > 1) {
        setSlice(nSlices);
        run("Delete Slice");
    }
    run("8-bit");
    setThreshold(1, 255);
    run("Convert to Mask");

    roiManager("reset");
    run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity add");
    close();

    n = roiManager("count");
    if (n == 0) {
        exit("No objects found in mask.\n"
            + "Check that nuclei are white on black background.");
    }

    // Measure on original image
    run("Clear Results");
    for (i = 0; i < n; i++) {
        selectImage(imgID);
        roiManager("select", i);
        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            addToResultsTable("Nucleus_" + (i + 1), results);
        }
    }
    roiManager("deselect");
    showCompletionMessage(n, "Binary Mask");
}

// =========================================================================
// MODE 4: Auto-Segmentation
// =========================================================================
function modeAutoSegmentation(imgID, excludeEdges, savedCh, savedSl, savedFr) {
    selectImage(imgID);

    // Extract only the current single 2D plane from any image type
    if (nSlices > 1) {
        run("Duplicate...", "title=_NM_temp_autoseg_ duplicate channels=" + savedCh
            + " slices=" + savedSl + " frames=" + savedFr);
    } else {
        run("Duplicate...", "title=_NM_temp_autoseg_");
    }
    tempID = getImageID();

    // Safety: ensure single 2D slice
    while (nSlices > 1) {
        setSlice(nSlices);
        run("Delete Slice");
    }

    // Ensure grayscale
    if (bitDepth() == 24) {
        run("8-bit");
    }

    // Preprocessing
    run("Gaussian Blur...", "sigma=2");

    // Threshold
    setAutoThreshold(AUTO_THRESH_METHOD + " dark");
    run("Convert to Mask");

    // Cleanup
    run("Fill Holes");
    run("Open");

    // Detect nuclei
    roiManager("reset");
    if (excludeEdges) {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity exclude add");
    } else {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity add");
    }
    close();  // close temp segmentation image

    n = roiManager("count");
    if (n == 0) {
        exit("No nuclei found.\n"
            + "Try a different threshold method or lower min area.");
    }

    // Return to original image at the SAME slice the user was viewing
    selectImage(imgID);
    if (nSlices > 1) {
        Stack.setPosition(savedCh, savedSl, savedFr);
    }
    roiManager("Show All");

    // Measure on original image at the saved slice position
    run("Clear Results");
    for (i = 0; i < n; i++) {
        selectImage(imgID);
        if (nSlices > 1) {
            Stack.setPosition(savedCh, savedSl, savedFr);
        }
        roiManager("select", i);
        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            addToResultsTable("Auto_" + (i + 1), results);
        }
    }
    roiManager("deselect");

    // Final: restore original view position
    selectImage(imgID);
    if (nSlices > 1) {
        Stack.setPosition(savedCh, savedSl, savedFr);
    }

    showCompletionMessage(n, "Auto-Segmentation");
}

// =========================================================================
// CORE: Extract pixel intensities from current selection
// =========================================================================
function getPixelsFromSelection(imgID) {
    selectImage(imgID);
    getSelectionBounds(rx, ry, rw, rh);

    pixelValues = newArray(0);
    for (yy = ry; yy < ry + rh; yy++) {
        for (xx = rx; xx < rx + rw; xx++) {
            if (selectionContains(xx, yy)) {
                v = getPixel(xx, yy);
                pixelValues = Array.concat(pixelValues, v);
            }
        }
    }
    return pixelValues;
}

// =========================================================================
// CORE: Compute CV, 1-Gini, DSI
//
// Formulas:
//   CV  = sigma / mu                          (on raw pixel intensities)
//   Normalized: x_i = (pixel_i - min) / (max - min)
//   DSI = count(x_i > tau) / N                (fraction above threshold)
//   Gini = [2 * sum(i * x_sorted_i)] / [N * sum(x)] - (N+1)/N
//   1-Gini = 1 - Gini
// =========================================================================
function computeMetrics(pixels) {
    n = pixels.length;

    // --- CV on raw intensities ---
    Array.getStatistics(pixels, minVal, maxVal, mu, sigma);
    if (mu != 0) {
        cv = sigma / mu;
    } else {
        cv = 0;
    }

    // --- Min-max normalization ---
    range = maxVal - minVal;
    normPixels = newArray(n);
    if (range > 0) {
        for (i = 0; i < n; i++) {
            normPixels[i] = (pixels[i] - minVal) / range;
        }
    }

    // --- DSI ---
    countAbove = 0;
    for (i = 0; i < n; i++) {
        if (normPixels[i] > DSI_THRESHOLD) {
            countAbove++;
        }
    }
    dsi = countAbove / n;

    // --- Gini coefficient ---
    sortedNorm = Array.copy(normPixels);
    Array.sort(sortedNorm);

    sumAll = 0;
    weightedSum = 0;
    for (i = 0; i < n; i++) {
        sumAll += sortedNorm[i];
        weightedSum += (i + 1) * sortedNorm[i];
    }

    gini = 0;
    if (sumAll > 0) {
        gini = ((2 * weightedSum) / (n * sumAll)) - ((n + 1) / n);
    }
    oneMinusGini = 1 - gini;

    // Return: [CV, 1-Gini, DSI, mean, stddev, n_pixels]
    return newArray(cv, oneMinusGini, dsi, mu, sigma, n);
}

// =========================================================================
// OUTPUT: Log window
// =========================================================================
function printResults(label, r) {
    print("\\n============ NucMetrics ============");
    print("Region:       " + label);
    print("N pixels:     " + d2s(r[5], 0));
    print("Mean (raw):   " + d2s(r[3], 2));
    print("StdDev (raw): " + d2s(r[4], 2));
    print("------------------------------------");
    print("CV:           " + d2s(r[0], 4));
    print("1-Gini:       " + d2s(r[1], 4));
    print("DSI (tau=" + d2s(DSI_THRESHOLD, 2) + "): " + d2s(r[2], 4));
    print("====================================");
}

// =========================================================================
// OUTPUT: Results Table
// =========================================================================
function addToResultsTable(label, r) {
    row = nResults;
    setResult("Label", row, label);
    setResult("N_pixels", row, d2s(r[5], 0));
    setResult("Mean_raw", row, d2s(r[3], 2));
    setResult("StdDev_raw", row, d2s(r[4], 2));
    setResult("CV", row, d2s(r[0], 4));
    setResult("1-Gini", row, d2s(r[1], 4));
    setResult("DSI", row, d2s(r[2], 4));
    updateResults();
}

// =========================================================================
// OUTPUT: Completion message
// =========================================================================
function showCompletionMessage(n, modeName) {
    print("\\n[NucMetrics] Processed " + n + " nuclei (" + modeName + ").");
    print("Results are in the Results table.");
    print("Cite: Kang et al. (2026) bioRxiv doi:10.64898/2026.03.30.715467");
}
