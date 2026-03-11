import os
import yaml
import numpy as np
import pandas as pd
from skimage import io, measure

# ---------------------------------------------------------
# Constants
# ---------------------------------------------------------
PIXEL_SIZE_UM = 0.122249  # micrometers per pixel
AREA_SCALE = PIXEL_SIZE_UM ** 2

# ---------------------------------------------------------
# Load Paths
# ---------------------------------------------------------
CONFIG_FILE = "local_config.yml"
with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

PARENT_DIR = "/Users/lukeishag/Documents/PERL/CAP2/manually_curated"
#os.path.join(local_config["raw_data_path"],"../60x/derivatives/segmentations/manually_curated")

OUTPUT_CSV = os.path.join(
    local_config["raw_data_path"],
    "..", "60x", "derivatives", "summary", "fascicle_data.csv"
)

# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

def load_mask(path):
    """Load a binary mask and ensure it's boolean."""
    img = io.imread(path)
    return img.astype(bool)

def extract_region_metrics(mask):
    """
    Extract geometric metrics for a single binary mask.
    Converts all pixel-based metrics into micrometers.
    """
    labeled = measure.label(mask)
    props = measure.regionprops(labeled)

    if len(props) == 0:
        return None

    r = max(props, key=lambda x: x.area)

    area_um2 = r.area * AREA_SCALE
    perimeter_um = r.perimeter * PIXEL_SIZE_UM

    # area-derived diameter
    area_diam_um = 2 * np.sqrt(area_um2 / np.pi)

    return {
        "area_um2": area_um2,
        "perimeter_um": perimeter_um,
        "major_axis_um": r.major_axis_length * PIXEL_SIZE_UM,
        "minor_axis_um": r.minor_axis_length * PIXEL_SIZE_UM,
        "eccentricity": r.eccentricity,
        "orientation": r.orientation,
        "centroid_row_um": r.centroid[0] * PIXEL_SIZE_UM,
        "centroid_col_um": r.centroid[1] * PIXEL_SIZE_UM,
        "area_diam_um": area_diam_um,
    }

# ---------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------

rows = []

for root, dirs, files in os.walk(PARENT_DIR):

    # Identify fascicle masks
    i_files = sorted([f for f in files if f.startswith("i") and f.endswith(".tif") and f[1].isdigit()])
    if len(i_files) == 0:
        continue

    subject = os.path.basename(root)

    # -----------------------------------------------------
    # Load axon mask (f.tif)
    # -----------------------------------------------------
    f_path = os.path.join(root, "f.tif")
    axon_area_um2 = 0
    if os.path.exists(f_path):
        f_mask = load_mask(f_path)
        axon_area_um2 = np.sum(f_mask) * AREA_SCALE

    # -----------------------------------------------------
    # Load outer boundary mask (o.tif)
    # -----------------------------------------------------
    o_path = os.path.join(root, "o.tif")
    outer_metrics = None

    if os.path.exists(o_path):
        o_mask = load_mask(o_path)
        outer_metrics = extract_region_metrics(o_mask)

        if outer_metrics is not None:
            rows.append({
                "subject": subject,
                "mask_type": "o",
                **outer_metrics
            })

    # -----------------------------------------------------
    # Process each fascicle mask i1.tif, i2.tif, ...
    # -----------------------------------------------------
    fascicle_areas_um2 = []
    fascicle_masks = []

    for fname in i_files:
        i_path = os.path.join(root, fname)
        mask_i = load_mask(i_path)
        fascicle_masks.append(mask_i)

        metrics = extract_region_metrics(mask_i)
        if metrics is None:
            continue

        fascicle_areas_um2.append(metrics["area_um2"])

        rows.append({
            "subject": subject,
            "mask_type": "i",
            "fascicle_id": fname.replace(".tif", ""),
            **metrics
        })

    # -----------------------------------------------------
    # Whole-nerve ellipse fit (combine all masks)
    # -----------------------------------------------------
    combined_mask = np.zeros_like(fascicle_masks[0], dtype=bool)

    for m in fascicle_masks:
        combined_mask |= m

    if outer_metrics is not None:
        combined_mask |= load_mask(o_path)

    whole_metrics = extract_region_metrics(combined_mask)

    if whole_metrics is not None:
        rows.append({
            "subject": subject,
            "mask_type": "whole_nerve",
            **whole_metrics
        })

    # -----------------------------------------------------
    # Compute total area of all fascicles
    # -----------------------------------------------------
    total_area_um2 = sum(fascicle_areas_um2)

    # -----------------------------------------------------
    # Compute axon-to-fascicle area ratio
    # -----------------------------------------------------
    area_ratio = axon_area_um2 / total_area_um2 if total_area_um2 > 0 else np.nan

    # -----------------------------------------------------
    # Add summary row for the subject
    # -----------------------------------------------------
    summary_row = {
        "subject": subject,
        "mask_type": "summary",
        "num_fascicles": len(fascicle_areas_um2),
        "total_area_um2": total_area_um2,
        "axon_area_um2": axon_area_um2,
        "area_ratio": area_ratio,
    }

    if outer_metrics is not None:
        summary_row.update({
            "outer_area_um2": outer_metrics["area_um2"],
            "outer_perimeter_um": outer_metrics["perimeter_um"],
            "outer_thickness_um": outer_metrics["area_um2"] / outer_metrics["perimeter_um"]
                if outer_metrics["perimeter_um"] > 0 else np.nan,
            "outer_area_diam_um": outer_metrics["area_diam_um"]
        })

    if whole_metrics is not None:
        summary_row.update({
            "whole_area_um2": whole_metrics["area_um2"],
            "whole_area_diam_um": whole_metrics["area_diam_um"],
            "whole_major_axis_um": whole_metrics["major_axis_um"],
            "whole_minor_axis_um": whole_metrics["minor_axis_um"],
            "whole_eccentricity": whole_metrics["eccentricity"]
        })

    rows.append(summary_row)

# ---------------------------------------------------------
# SAVE CSV
# ---------------------------------------------------------

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)

print(f"Saved fascicle metrics for {len(df)} rows to {OUTPUT_CSV}")
