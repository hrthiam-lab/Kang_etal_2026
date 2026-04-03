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
//   3. Binary Mask        - use an open binary mask image/stack to define nuclei
//   4. Auto-Generate Binary Mask - threshold + cleanup + mask review
//                           (current slice or entire stack)
//
// Installation:
//   1. Copy this file to: Fiji.app/macros/toolsets/
//      (or ImageJ/macros/toolsets/)
//   2. Restart Fiji/ImageJ
//   3. Click ">>" on the toolbar, select "NucMetrics_Toolset"
//   4. Click the NucMetrics icon to run
//
// If you use this tool, please cite:
//   Kang M, Cabral AT, Sawant M, Thiam HR (2026).
//   "Benchmarking three simple DNA staining-based image metrics for
//   live-cell tracking of chromatin organization."
//   bioRxiv. https://doi.org/10.64898/2026.03.30.715467
//
// Recent updates:
//   - Current Selection mode: guides user to draw ROI if none exists
//   - Binary Mask mode: current-slice and whole-stack support for matching mask stacks
//   - Auto-Generate Binary Mask: generate, review, then optionally compute
//   - Hyperstack auto-mask: choose current C/Z across T or current C/T across Z
//   - Improved ImageJ (non-Fiji) compatibility for stack handling
//
// Created: 2026-04-01 | Updated: 2026-04-03
// Authors: Minwoo Kang Ph.D.
//          HR Thiam Lab, Stanford University
// License: MIT
// =========================================================================

// ---- Global parameters ----
var DSI_THRESHOLD = 0.3;
var MIN_AREA = 30;
var AUTO_THRESH_METHOD = "Li";

// =========================================================================
// TOOLBAR ICON
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

macro "NucMetrics [n1]" {
    runNucMetrics();
}

// =========================================================================
// HELPER: Get current slice position (ImageJ + Fiji compatible)
// Stack.getPosition() is Fiji-only; this falls back to getSliceNumber()
// =========================================================================
function getCurrentStackPosition() {
    // Returns array [channel, slice, frame]
    ch = 1; sl = 1; fr = 1;
    getDimensions(w, h, channels, slices, frames);
    if (is("hyperstack")) {
        Stack.getPosition(ch, sl, fr);
    } else if (nSlices > 1) {
        // Plain stack: current plane is represented by the slice index
        sl = getSliceNumber();
    }
    return newArray(ch, sl, fr);
}

// =========================================================================
// HELPER: Set stack position (ImageJ + Fiji compatible)
// =========================================================================
function setStackPosition(ch, sl, fr) {
    getDimensions(w, h, channels, slices, frames);
    if (is("hyperstack")) {
        Stack.setPosition(ch, sl, fr);
    } else if (nSlices > 1) {
        setSlice(sl);
    }
}

// =========================================================================
// HELPER: Number of planes/time points to iterate
// For hyperstacks, iterate over frames at the current C/Z.
// For plain stacks, iterate over slices.
// =========================================================================
function getSeriesLength(imgID) {
    selectImage(imgID);
    getDimensions(w, h, channels, slices, frames);
    if (is("hyperstack")) {
        return frames;
    } else if (nSlices > 1) {
        return nSlices;
    } else {
        return 1;
    }
}

// =========================================================================
// HELPER: Duplicate single plane (ImageJ + Fiji compatible)
// =========================================================================
function duplicateSinglePlane(title, ch, sl, fr) {
    getDimensions(w, h, channels, slices, frames);

    if (is("hyperstack")) {
        cRange = ch + "-" + ch;
        zRange = sl + "-" + sl;
        tRange = fr + "-" + fr;
        run("Duplicate...", "title=" + title + " duplicate channels=" + cRange
            + " slices=" + zRange + " frames=" + tRange);
    } else if (nSlices > 1) {
        // Plain stack: explicitly duplicate only the requested plane.
        // Using range=sl-sl avoids accidentally duplicating the whole stack.
        setSlice(sl);
        run("Duplicate...", "title=" + title + " duplicate range=" + sl + "-" + sl);
    } else {
        run("Duplicate...", "title=" + title);
    }

    // Safety: ensure the duplicate is a single 2D image.
    while (nSlices > 1) {
        setSlice(nSlices);
        run("Delete Slice");
    }
}

// =========================================================================
// MAIN
// =========================================================================
function runNucMetrics() {

    if (nImages == 0) {
        exit("No image open.\nPlease open a DNA-stained image first.");
    }

    originalID = getImageID();

    // Save current position (compatible with ImageJ and Fiji)
    pos = getCurrentStackPosition();
    origChannel = pos[0];
    origSlice = pos[1];
    origFrame = pos[2];

    // --- STEP 1: Mode + DSI threshold ---
    modes = newArray(
        "Current Selection",
        "ROI Manager",
        "Binary Mask",
        "Auto-Generate Binary Mask"
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

    // --- STEP 2: Mode-specific dialog and execution ---
    if (mode == "Auto-Generate Binary Mask") {
        getDimensions(w, h, channels, slices, frames);
        imgIsHyper = is("hyperstack");
        isStack = false;
        autoSeriesMode = "single";
        autoSeriesLabel = "S";
        showAxisChoice = false;

        if (imgIsHyper) {
            if (frames > 1 || slices > 1) {
                isStack = true;
            }
            if (frames > 1 && slices > 1) {
                showAxisChoice = true;
            } else if (frames > 1) {
                autoSeriesMode = "frames";
                autoSeriesLabel = "T";
            } else if (slices > 1) {
                autoSeriesMode = "slices";
                autoSeriesLabel = "Z";
            }
        } else if (nSlices > 1) {
            isStack = true;
            autoSeriesMode = "stack";
            autoSeriesLabel = "S";
        }

        Dialog.create("NucMetrics - Auto-Generate Binary Mask");
        Dialog.addMessage("This mode now generates a binary mask first.\nIt shows the mask in a new window, and computes metrics only after you confirm the mask.");
        Dialog.addChoice("Threshold method:",
            newArray("Li", "Otsu", "Triangle"), AUTO_THRESH_METHOD);
        Dialog.addNumber("Min nucleus area (px):", MIN_AREA);
        Dialog.addCheckbox("Exclude nuclei touching edges", true);
        if (channels > 1) {
            Dialog.addMessage("Active channel only: C" + origChannel);
        }
        if (isStack) {
            Dialog.addMessage("--- Stack options ---");
            Dialog.addCheckbox("Generate mask for an entire series", false);
            if (showAxisChoice) {
                Dialog.addChoice("Series axis:",
                    newArray("Current channel + current Z across time (T)",
                             "Current channel + current time across Z"),
                    "Current channel + current Z across time (T)");
            } else if (imgIsHyper && frames > 1) {
                Dialog.addMessage("Series axis: current channel + current Z across time (T)");
            } else if (imgIsHyper && slices > 1) {
                Dialog.addMessage("Series axis: current channel + current time across Z");
            } else {
                Dialog.addMessage("Series axis: plain stack order (S)");
            }
        }
        Dialog.show();
        AUTO_THRESH_METHOD = Dialog.getChoice();
        MIN_AREA = Dialog.getNumber();
        excludeEdges = Dialog.getCheckbox();
        processAllSlices = false;
        if (isStack) {
            processAllSlices = Dialog.getCheckbox();
            if (showAxisChoice) {
                axisChoice = Dialog.getChoice();
                if (axisChoice == "Current channel + current Z across time (T)") {
                    autoSeriesMode = "frames";
                    autoSeriesLabel = "T";
                } else {
                    autoSeriesMode = "slices";
                    autoSeriesLabel = "Z";
                }
            }
        }

        modeAutoSegGenerateMaskAndMaybeCompute(originalID, processAllSlices,
            excludeEdges, origChannel, origSlice, origFrame, autoSeriesMode, autoSeriesLabel);

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
        Dialog.addMessage("The binary mask image or mask stack must already be open.\n"
            + "Select one of the currently open mask windows below.");
        Dialog.addChoice("Mask image/stack:", maskList, maskList[0]);
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

        processMaskAll = false;
        imgSeriesLength = getSeriesLength(originalID);
        maskSeriesLength = getSeriesLength(maskID);
        if (imgSeriesLength > 1 || maskSeriesLength > 1) {
            Dialog.create("NucMetrics - Binary Mask Stack Options");
            Dialog.addMessage("Current image planes/time points: " + imgSeriesLength + "\n"
                + "Mask planes/time points: " + maskSeriesLength + "\n\n"
                + "Whole-stack binary-mask mode measures each corresponding\n"
                + "frame/slice using the provided mask stack.");
            Dialog.addCheckbox("Process entire stack (matching planes/time points)", false);
            Dialog.show();
            processMaskAll = Dialog.getCheckbox();
        }

        if (processMaskAll) {
            modeBinaryMaskAllSlices(originalID, maskID, origChannel, origSlice, origFrame);
        } else {
            modeBinaryMask(originalID, maskID);
        }

    } else if (mode == "ROI Manager") {
        modeROIManager(originalID);

    } else {
        modeCurrentSelection(originalID);
    }
}

// =========================================================================
// MODE 1: Current Selection (with ROI guidance)
// =========================================================================
function modeCurrentSelection(imgID) {
    selectImage(imgID);

    // If no ROI exists, prompt user to draw one
    if (selectionType() == -1) {
        // Set freehand tool
        setTool("freehand");
        waitForUser("NucMetrics - Draw ROI",
            "No selection found.\n \n"
            + "Please draw an ROI around a nucleus using\n"
            + "the Freehand tool (now active), then click OK.\n \n"
            + "Tip: You can also use Polygon, Oval, or\n"
            + "any other selection tool from the toolbar.");
    }

    // Check again after user interaction
    if (selectionType() == -1) {
        exit("No selection drawn. Operation cancelled.");
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
        exit("ROI Manager is empty or not open.\nOpen ROI Manager, add nuclear ROIs, then try again.");
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

    selectImage(maskID);
    pos = getCurrentStackPosition();
    duplicateSinglePlane("_NM_temp_mask_", pos[0], pos[1], pos[2]);
    tempMaskID = getImageID();
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
// MODE 3b: Binary Mask - ENTIRE STACK
// =========================================================================
function modeBinaryMaskAllSlices(imgID, maskID, savedCh, savedSl, savedFr) {
    // Check XY dimensions
    selectImage(imgID);
    w1 = getWidth(); h1 = getHeight();
    getDimensions(iw, ih, ichannels, islices, iframes);
    imgIsHyper = is("hyperstack");

    selectImage(maskID);
    w2 = getWidth(); h2 = getHeight();
    getDimensions(mw, mh, mchannels, mslices, mframes);
    maskIsHyper = is("hyperstack");

    if (w1 != w2 || h1 != h2) {
        exit("Dimension mismatch.\n"
            + "Image: " + w1 + "x" + h1 + "\n"
            + "Mask:  " + w2 + "x" + h2);
    }

    // Decide iteration axis for image
    imgMode = "stack";
    imgLen = 1;
    axisLabel = "S";
    if (imgIsHyper) {
        if (iframes > 1) {
            imgMode = "frames";
            imgLen = iframes;
            axisLabel = "T";
        } else if (islices > 1) {
            imgMode = "slices";
            imgLen = islices;
            axisLabel = "Z";
        } else {
            imgLen = 1;
        }
    } else if (nSlices > 1) {
        imgMode = "stack";
        imgLen = nSlices;
        axisLabel = "S";
    }

    // Decide iteration axis for mask
    selectImage(maskID);
    maskMode = "stack";
    maskLen = 1;
    if (maskIsHyper) {
        if (mframes > 1) {
            maskMode = "frames";
            maskLen = mframes;
        } else if (mslices > 1) {
            maskMode = "slices";
            maskLen = mslices;
        } else {
            maskLen = 1;
        }
    } else if (nSlices > 1) {
        maskMode = "stack";
        maskLen = nSlices;
    }

    if (imgLen != maskLen) {
        exit("Image and mask stack lengths do not match.\n"
            + "Image length: " + imgLen + "\n"
            + "Mask length:  " + maskLen);
    }

    if (imgLen <= 1) {
        modeBinaryMask(imgID, maskID);
        return;
    }

    run("Clear Results");
    setBatchMode(true);
    measuredPlanes = 0;
    skippedPlanes = 0;
    roiManager("reset");

    for (p = 1; p <= imgLen; p++) {
        showProgress(p, imgLen);
        roiManager("reset");

        // Get selection from the matching mask plane
        selectImage(maskID);
        if (maskMode == "frames") {
            setStackPosition(savedCh, savedSl, p);
            duplicateSinglePlane("_NM_temp_mask_stack_", savedCh, savedSl, p);
        } else if (maskMode == "slices") {
            setStackPosition(savedCh, p, savedFr);
            duplicateSinglePlane("_NM_temp_mask_stack_", savedCh, p, savedFr);
        } else {
            setSlice(p);
            duplicateSinglePlane("_NM_temp_mask_stack_", 1, p, 1);
        }

        run("8-bit");
        setThreshold(1, 255);
        run("Convert to Mask");
        run("Create Selection");

        if (selectionType() == -1) {
            skippedPlanes++;
            close();
            continue;
        }
        roiManager("Add");
        close();

        // Measure on the matching image plane using the ROI transferred from the mask plane
        selectImage(imgID);
        if (imgMode == "frames") {
            setStackPosition(savedCh, savedSl, p);
        } else if (imgMode == "slices") {
            setStackPosition(savedCh, p, savedFr);
        } else {
            setSlice(p);
        }
        roiManager("Select", 0);

        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            label = axisLabel + IJ.pad(p, 3) + "_N1";
            addToResultsTable(label, results);
            measuredPlanes++;
        } else {
            skippedPlanes++;
        }
    }
    roiManager("reset");

    setBatchMode(false);
    selectImage(imgID);
    setStackPosition(savedCh, savedSl, savedFr);

    print("\n[NucMetrics] Processed " + imgLen + " planes using a binary mask stack.");
    print("Measured planes: " + measuredPlanes + "; skipped planes: " + skippedPlanes + ".");
    print("Label format: " + axisLabel + "###_N1 (plane index + nucleus 1).");
    print("Results are in the Results table.");
    print("Cite: Kang et al. (2026) bioRxiv doi:10.64898/2026.03.30.715467");
}

// =========================================================================
// MODE 3 helper: Binary Mask whole-series measurement with a forced axis
// Used by Option 4 after generating a reviewable mask stack.
// =========================================================================
function modeBinaryMaskAllSlicesForced(imgID, maskID, savedCh, savedSl, savedFr, imgMode, axisLabel) {
    selectImage(imgID);
    w1 = getWidth(); h1 = getHeight();
    getDimensions(iw, ih, ichannels, islices, iframes);

    selectImage(maskID);
    w2 = getWidth(); h2 = getHeight();
    maskLen = getSeriesLength(maskID);

    if (w1 != w2 || h1 != h2) {
        exit("Dimension mismatch.\nImage: " + w1 + "x" + h1 + "\nMask:  " + w2 + "x" + h2);
    }

    imgLen = 1;
    if (imgMode == "frames") {
        imgLen = iframes;
    } else if (imgMode == "slices") {
        imgLen = islices;
    } else if (imgMode == "stack") {
        imgLen = nSlices;
    }

    if (imgLen != maskLen) {
        exit("Image and generated mask lengths do not match.\nImage length: " + imgLen + "\nMask length:  " + maskLen);
    }

    if (imgLen <= 1) {
        modeBinaryMask(imgID, maskID);
        return;
    }

    run("Clear Results");
    setBatchMode(true);
    measuredPlanes = 0;
    skippedPlanes = 0;
    roiManager("reset");

    for (p = 1; p <= imgLen; p++) {
        showProgress(p, imgLen);
        roiManager("reset");

        selectImage(maskID);
        setSlice(p);
        duplicateSinglePlane("_NM_temp_mask_stack_", 1, p, 1);
        run("8-bit");
        setThreshold(1, 255);
        run("Convert to Mask");
        run("Create Selection");

        if (selectionType() == -1) {
            skippedPlanes++;
            close();
            continue;
        }
        roiManager("Add");
        close();

        selectImage(imgID);
        if (imgMode == "frames") {
            setStackPosition(savedCh, savedSl, p);
        } else if (imgMode == "slices") {
            setStackPosition(savedCh, p, savedFr);
        } else {
            setSlice(p);
        }
        roiManager("Select", 0);

        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            label = axisLabel + IJ.pad(p, 3) + "_N1";
            addToResultsTable(label, results);
            measuredPlanes++;
        } else {
            skippedPlanes++;
        }
    }
    roiManager("reset");

    setBatchMode(false);
    selectImage(imgID);
    setStackPosition(savedCh, savedSl, savedFr);

    print("[NucMetrics] Processed " + imgLen + " planes using a generated binary mask series.");
    print("Series axis: " + axisLabel + " (current channel fixed at C" + savedCh + ").");
    print("Measured planes: " + measuredPlanes + "; skipped planes: " + skippedPlanes + ".");
    print("Label format: " + axisLabel + "###_N1 (plane index + nucleus 1).");
    print("Results are in the Results table.");
    print("Cite: Kang et al. (2026) bioRxiv doi:10.64898/2026.03.30.715467");
}

// =========================================================================
// MODE 4: Auto-Generate Binary Mask -> generate binary mask, review, then compute
// =========================================================================
function modeAutoSegGenerateMaskAndMaybeCompute(imgID, processAllSlices, excludeEdges,
    savedCh, savedSl, savedFr, seriesMode, seriesLabel) {

    if (processAllSlices) {
        maskID = generateAutoMaskStack(imgID, excludeEdges, savedCh, savedSl, savedFr, seriesMode);
    } else {
        maskID = generateAutoMaskSingle(imgID, excludeEdges, savedCh, savedSl, savedFr);
    }

    selectImage(maskID);
    if (processAllSlices) {
        print("\n[NucMetrics] Generated binary mask stack: " + getTitle());
    } else {
        print("\n[NucMetrics] Generated binary mask image: " + getTitle());
    }

    showMessageWithCancel("NucMetrics - Review Auto Mask",
        "A binary mask has been generated in a new window.\n\n"
        + "Review the mask and click OK to compute NucMetrics using this mask.\n"
        + "Click Cancel to stop and keep the generated mask only.");

    if (processAllSlices) {
        modeBinaryMaskAllSlicesForced(imgID, maskID, savedCh, savedSl, savedFr, seriesMode, seriesLabel);
    } else {
        modeBinaryMask(imgID, maskID);
    }
}

// =========================================================================
// HELPER: Create a single-slice binary mask from the current image plane
// =========================================================================
function generateAutoMaskSingle(imgID, excludeEdges, savedCh, savedSl, savedFr) {
    selectImage(imgID);
    baseTitle = getTitle();
    duplicateSinglePlane("_NM_autoMaskTemp_", savedCh, savedSl, savedFr);
    tempID = getImageID();
    processCurrentTempToBinaryMask(excludeEdges);
    rename(buildAutoMaskTitle(baseTitle, false));
    return getImageID();
}

// =========================================================================
// HELPER: Create a plain binary mask stack matching the processed series length
// =========================================================================
function generateAutoMaskStack(imgID, excludeEdges, savedCh, savedSl, savedFr, seriesMode) {
    selectImage(imgID);
    baseTitle = getTitle();
    getDimensions(w, h, channels, slices, frames);

    processMode = seriesMode;
    totalPlanes = 1;
    if (is("hyperstack")) {
        if (processMode == "frames") {
            totalPlanes = frames;
        } else if (processMode == "slices") {
            totalPlanes = slices;
        } else {
            totalPlanes = 1;
        }
    } else if (nSlices > 1) {
        processMode = "stack";
        totalPlanes = nSlices;
    }

    if (totalPlanes <= 1) {
        return generateAutoMaskSingle(imgID, excludeEdges, savedCh, savedSl, savedFr);
    }

    // Build a reduced working series first.
    // This is more reliable for multi-channel / multi-Z hyperstacks than
    // extracting one plane at a time from the original CxZxT image.
    workTitle = "_NM_workSeries_";
    selectImage(imgID);
    if (is("hyperstack")) {
        if (processMode == "frames") {
            run("Duplicate...", "title=" + workTitle + " duplicate channels=" + savedCh + "-" + savedCh
                + " slices=" + savedSl + "-" + savedSl + " frames=1-" + frames);
        } else if (processMode == "slices") {
            run("Duplicate...", "title=" + workTitle + " duplicate channels=" + savedCh + "-" + savedCh
                + " slices=1-" + slices + " frames=" + savedFr + "-" + savedFr);
        }
        if (is("hyperstack")) {
            run("Hyperstack to Stack");
        }
    } else {
        run("Duplicate...", "title=" + workTitle + " duplicate");
    }
    workID = getImageID();

    // Refresh the actual number of planes in the reduced series.
    selectImage(workID);
    workLen = nSlices;
    if (workLen <= 1) {
        close();
        return generateAutoMaskSingle(imgID, excludeEdges, savedCh, savedSl, savedFr);
    }

    outTitle = buildAutoMaskTitle(baseTitle, true);
    newImage(outTitle, "8-bit black", w, h, workLen);
    outID = getImageID();

    setBatchMode(true);
    for (p = 1; p <= workLen; p++) {
        showProgress(p, workLen);
        selectImage(workID);
        setSlice(p);
        duplicateSinglePlane("_NM_autoMaskTemp_", 1, p, 1);

        tempID = getImageID();
        processCurrentTempToBinaryMask(excludeEdges);
        run("Select All");
        run("Copy");

        selectImage(outID);
        setSlice(p);
        run("Paste");
        run("Select None");

        selectImage(tempID);
        close();
    }
    setBatchMode(false);

    selectImage(workID);
    close();

    selectImage(outID);
    setSlice(1);
    return outID;
}

// =========================================================================
// HELPER: Segment the current temporary image and leave it as a binary mask
// =========================================================================
function processCurrentTempToBinaryMask(excludeEdges) {
    if (bitDepth() == 24) run("8-bit");

    run("Gaussian Blur...", "sigma=2");
    setAutoThreshold(AUTO_THRESH_METHOD + " dark");
    run("Convert to Mask");
    run("Fill Holes");
    run("Open");

    roiManager("reset");
    if (excludeEdges) {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity exclude add");
    } else {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity add");
    }

    run("Select All");
    setForegroundColor(0, 0, 0);
    run("Fill");
    run("Select None");

    n = roiManager("count");
    if (n > 0) {
        setForegroundColor(255, 255, 255);
        for (i = 0; i < n; i++) {
            roiManager("select", i);
            run("Fill");
        }
    }
    roiManager("reset");
    run("Select None");
}

// =========================================================================
// HELPER: Build a predictable title for generated auto-mask images
// =========================================================================
function buildAutoMaskTitle(baseTitle, isStackMask) {
    if (isStackMask) {
        return baseTitle + "_AutoMaskStack";
    } else {
        return baseTitle + "_AutoMask";
    }
}

// =========================================================================
// MODE 4a: Auto-Generate Binary Mask - SINGLE SLICE
// =========================================================================
function modeAutoSegSingleSlice(imgID, excludeEdges, savedCh, savedSl, savedFr) {
    selectImage(imgID);
    duplicateSinglePlane("_NM_temp_autoseg_", savedCh, savedSl, savedFr);
    tempID = getImageID();

    if (bitDepth() == 24) {
        run("8-bit");
    }

    // Segment
    run("Gaussian Blur...", "sigma=2");
    setAutoThreshold(AUTO_THRESH_METHOD + " dark");
    run("Convert to Mask");
    run("Fill Holes");
    run("Open");

    roiManager("reset");
    if (excludeEdges) {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity exclude add");
    } else {
        run("Analyze Particles...", "size=" + MIN_AREA + "-Infinity add");
    }
    close();

    n = roiManager("count");
    if (n == 0) {
        exit("No nuclei found.\n"
            + "Try a different threshold method or lower min area.");
    }

    // Return to original and restore position
    selectImage(imgID);
    setStackPosition(savedCh, savedSl, savedFr);
    roiManager("Show All");

    // Measure
    run("Clear Results");
    for (i = 0; i < n; i++) {
        selectImage(imgID);
        setStackPosition(savedCh, savedSl, savedFr);
        roiManager("select", i);
        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            addToResultsTable("Auto_" + (i + 1), results);
        }
    }
    roiManager("deselect");
    selectImage(imgID);
    setStackPosition(savedCh, savedSl, savedFr);
    showCompletionMessage(n, "Auto-Generate Binary Mask (single slice)");
}

// =========================================================================
// MODE 4b: Auto-Generate Binary Mask - ENTIRE STACK
// =========================================================================
function modeAutoSegAllSlices(imgID, excludeEdges, savedCh, savedSl, savedFr) {
    selectImage(imgID);
    getDimensions(w, h, channels, slices, frames);

    // Decide which axis to iterate.
    processMode = "stack";
    totalPlanes = 1;
    axisLabel = "S";

    if (is("hyperstack")) {
        if (frames > 1) {
            processMode = "frames";
            totalPlanes = frames;
            axisLabel = "T";
        } else if (slices > 1) {
            processMode = "slices";
            totalPlanes = slices;
            axisLabel = "Z";
        } else {
            exit("This image does not contain multiple frames or slices.");
        }
    } else if (nSlices > 1) {
        processMode = "stack";
        totalPlanes = nSlices;
        axisLabel = "S";
    } else {
        exit("This image does not contain multiple frames or slices.");
    }

    // Whole-stack mode is intentionally limited to single-nucleus stacks.
    // It measures one auto-segmented nucleus per plane and does not perform tracking.
    run("Clear Results");
    roiManager("reset");
    setBatchMode(true);
    measuredPlanes = 0;
    skippedPlanes = 0;
    edgeWarningShown = false;

    for (p = 1; p <= totalPlanes; p++) {
        selectImage(imgID);
        showProgress(p, totalPlanes);

        if (processMode == "frames") {
            setStackPosition(savedCh, savedSl, p);
            duplicateSinglePlane("_NM_temp_stack_", savedCh, savedSl, p);
        } else if (processMode == "slices") {
            setStackPosition(savedCh, p, savedFr);
            duplicateSinglePlane("_NM_temp_stack_", savedCh, p, savedFr);
        } else {
            setSlice(p);
            duplicateSinglePlane("_NM_temp_stack_", 1, p, 1);
        }

        if (bitDepth() == 24) {
            run("8-bit");
        }

        // Segment the temporary single-plane image.
        run("Gaussian Blur...", "sigma=2");
        setAutoThreshold(AUTO_THRESH_METHOD + " dark");
        run("Convert to Mask");
        run("Fill Holes");
        run("Open");

        // Create one selection from the segmented mask.
        run("Create Selection");

        if (selectionType() == -1) {
            skippedPlanes++;
            close();
            continue;
        }

        getSelectionBounds(rx, ry, rw, rh);
        touchesEdge = (rx <= 0 || ry <= 0 || (rx + rw) >= w || (ry + rh) >= h);
        if (excludeEdges && touchesEdge) {
            if (!edgeWarningShown) {
                showMessage("NucMetrics Warning",
                    "At least one segmented nucleus touched the image border and was skipped.\n\n"
                    + "Whole-stack auto-segmentation is intended for single, fully visible nuclei.");
                edgeWarningShown = true;
            }
            skippedPlanes++;
            close();
            continue;
        }

        roiManager("reset");
        roiManager("add");
        close();

        // Measure on the original image at the matching plane.
        selectImage(imgID);
        if (processMode == "frames") {
            setStackPosition(savedCh, savedSl, p);
        } else if (processMode == "slices") {
            setStackPosition(savedCh, p, savedFr);
        } else {
            setSlice(p);
        }

        roiManager("select", 0);
        pixels = getPixelsFromSelection(imgID);
        if (pixels.length > 0) {
            results = computeMetrics(pixels);
            label = axisLabel + IJ.pad(p, 3) + "_N1";
            addToResultsTable(label, results);
            measuredPlanes++;
        } else {
            skippedPlanes++;
        }
    }

    setBatchMode(false);

    // Restore original position
    selectImage(imgID);
    setStackPosition(savedCh, savedSl, savedFr);
    roiManager("deselect");

    print("\n[NucMetrics] Processed " + totalPlanes + " planes in single-nucleus whole-stack mode.");
    print("Measured planes: " + measuredPlanes + "; skipped planes: " + skippedPlanes + ".");
    print("Label format: " + axisLabel + "###_N1 (plane index + nucleus 1).");
    print("Results are in the Results table.");
    print("Cite: Kang et al. (2026) bioRxiv doi:10.64898/2026.03.30.715467");
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
// =========================================================================
function computeMetrics(pixels) {
    n = pixels.length;

    Array.getStatistics(pixels, minVal, maxVal, mu, sigma);
    if (mu != 0) {
        cv = sigma / mu;
    } else {
        cv = 0;
    }

    range = maxVal - minVal;
    normPixels = newArray(n);
    for (i = 0; i < n; i++) {
        normPixels[i] = 0;
    }
    if (range > 0) {
        for (i = 0; i < n; i++) {
            normPixels[i] = (pixels[i] - minVal) / range;
        }
    }

    countAbove = 0;
    for (i = 0; i < n; i++) {
        if (normPixels[i] > DSI_THRESHOLD) {
            countAbove++;
        }
    }
    dsi = countAbove / n;

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
