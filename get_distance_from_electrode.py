import os
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from tifffile import imread  # SAFE TIFF LOADING

# Pixel size (µm per pixel)
PIXEL_SIZE_UM = 0.122249

# ---------------------------------------------------------
# Cuff diameter dictionaries
# ---------------------------------------------------------
stim_cuff_map = {
    "252 D05": 300,
    "262 D05": 300,
    "NCW-0.75-2-PT/IR-X(1-1)-1.5-00-300NE-Ban": 750,
    "NCW-1-2-PT/IR-X(1-1)-1.5-00-300-CONN": 1000,
    "NCX-1-2-500PR_1-1.5-300NE": 1000
}

# Histology Name → Expdate mapping
hist_to_expdate = {
    "39-1": "20220531", "39-2": "20220531",
    "40-1": "20220622",
    "41-1": "20220720", "41-2": "20220720",
    "43-1": "20221103",
    "44-1": "20221108", "44-2": "20221108",
    "45-1": "20221122",
    "46-1": "20221201", "46-2": "20221201",
    "47-1": "20230131",
    "48-1": "20230207", "48-2": "20230207"
}

# ---------------------------------------------------------
# Load config and CSVs
# ---------------------------------------------------------
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

RAW_DATA_PATH = local_config.get("raw_data_path")
BASE_SUMMARY_DIR = os.path.join(RAW_DATA_PATH, "..", "60x", "derivatives", "summary")

HIST_CSV = os.path.join(BASE_SUMMARY_DIR, "all_samples_fiber_diameters.csv")
FASCICLE_CSV = os.path.join(BASE_SUMMARY_DIR, "fascicle_data.csv")

hist_df = pd.read_csv(HIST_CSV)
fas_df = pd.read_csv(FASCICLE_CSV)

# ---------------------------------------------------------
# Convert all histology coordinates from µm → pixels
# ---------------------------------------------------------
hist_df["CentroidX_px"] = hist_df["CentroidX"] / PIXEL_SIZE_UM
hist_df["CentroidY_px"] = hist_df["CentroidY"] / PIXEL_SIZE_UM

fas_df["centroid_col_px"] = fas_df["centroid_col_um"] / PIXEL_SIZE_UM
fas_df["centroid_row_px"] = fas_df["centroid_row_um"] / PIXEL_SIZE_UM

# ---------------------------------------------------------
# Load metadata Excel to get stim_cuff
# ---------------------------------------------------------
metadata_path = str(os.path.join(RAW_DATA_PATH, "AnalysisScript_CAPs_Batch_parameters.xlsx"))
metadata_df = pd.read_excel(metadata_path)

# Normalize expdate_values to string for safe matching
metadata_df["expdate_values_str"] = metadata_df["expdate_values"].astype(str)

# Build expdate → stim cuff diameter lookup
expdate_to_stimdiam = {}
for _, mrow in metadata_df.iterrows():
    expdate_str = mrow["expdate_values_str"]
    stim_cuff = mrow["stim_cuff"]

    if stim_cuff not in stim_cuff_map:
        print(f"WARNING: stim_cuff '{stim_cuff}' not in stim_cuff_map, defaulting to 300 µm")

    expdate_to_stimdiam[expdate_str] = stim_cuff_map.get(stim_cuff, 300)

# ---------------------------------------------------------
# Add stim cuff diameter to hist_df
# ---------------------------------------------------------
diameters = []
for _, hrow in hist_df.iterrows():
    hist_name = hrow["Name"]
    expdate = hist_to_expdate.get(hist_name, None)

    if expdate is None:
        diameters.append(np.nan)
        continue

    diam_um = expdate_to_stimdiam.get(str(expdate), np.nan)
    diameters.append(diam_um)

hist_df["stim_cuff_diameter_um"] = diameters
hist_df["stim_cuff_radius_px"] = hist_df["stim_cuff_diameter_um"] / 2 / PIXEL_SIZE_UM

# ---------------------------------------------------------
# Build cache: subject → whole-nerve centroid (x, y) in pixels
# ---------------------------------------------------------
centroid_cache = {}

for subject in hist_df["Name"].unique():
    subset = fas_df[fas_df["subject"] == subject]

    whole_row = subset[subset["mask_type"] == "whole_nerve"]
    if whole_row.empty:
        print(f"Warning: No whole_nerve row found for subject {subject}")
        continue

    frow = whole_row.iloc[0]
    cx = frow["centroid_col_px"]
    cy = frow["centroid_row_px"]

    centroid_cache[subject] = (cx, cy)

# ---------------------------------------------------------
# Compute shortest distance from axon center to stim cuff circle (in µm)
# ---------------------------------------------------------
distance_to_stim_cuff = []

for idx, hrow in hist_df.iterrows():
    subject = hrow["Name"]

    # Skip if missing data
    if subject not in centroid_cache or np.isnan(hrow["stim_cuff_diameter_um"]):
        distance_to_stim_cuff.append(np.nan)
        continue

    # Nerve center in µm
    nerve_cx_px, nerve_cy_px = centroid_cache[subject]
    nerve_cx_um = nerve_cx_px * PIXEL_SIZE_UM
    nerve_cy_um = nerve_cy_px * PIXEL_SIZE_UM

    # Axon center in µm
    ax_x_um = hrow["CentroidX"]
    ax_y_um = hrow["CentroidY"]

    # Euclidean distance from axon to nerve center (µm)
    d_um = np.sqrt((ax_x_um - nerve_cx_um)**2 + (ax_y_um - nerve_cy_um)**2)

    # Stim cuff radius in µm
    cuff_radius_um = hrow["stim_cuff_diameter_um"] / 2.0

    # Shortest distance from axon to cuff boundary (µm)
    shortest_um = abs(d_um - cuff_radius_um)

    # Store ONLY this one value
    distance_to_stim_cuff.append(shortest_um)

# Only ONE new column is added to hist_df
hist_df["distance_to_stim_cuff_um"] = distance_to_stim_cuff

# ---------------------------------------------------------
# Save output
# ---------------------------------------------------------
LOCAL_OUTPUT_PATH = local_config["local_output_path"]
output_csv = os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv")

os.makedirs(LOCAL_OUTPUT_PATH, exist_ok=True)
hist_df.to_csv(output_csv, index=False)

print(f"Saved updated CSV to: {output_csv}")

# ---------------------------------------------------------
# Gut-check visualization
# ---------------------------------------------------------
subject_to_plot = "40-1"
hist_id = subject_to_plot
print("Hist ID:", hist_id)
print("Mapped expdate:", hist_to_expdate.get(hist_id))

# Build TIFF path automatically
IMAGE_PATH = f"/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/60x/raw/tiff/{subject_to_plot}.tif"

# Load TIFF safely
img = imread(IMAGE_PATH)

# Pull all histology rows for this subject
sub_hist = hist_df[hist_df["Name"] == subject_to_plot]

# Get nerve centroid (already in px)
nerve_cx_px, nerve_cy_px = centroid_cache[subject_to_plot]

# Get stim cuff radius for this subject (already in px)
subject_radius_px = sub_hist["stim_cuff_radius_px"].iloc[0]
print("stim_cuff_radius_px:", subject_radius_px)

# ---------------------------------------------------------
# Plot 1: centers + stim cuff circle
# ---------------------------------------------------------
plt.figure(figsize=(8, 8))
plt.imshow(img, cmap="gray")

# Plot axon centers
plt.scatter(
    sub_hist["CentroidX_px"],
    sub_hist["CentroidY_px"],
    s=5,
    c="red"
)

# Draw stim cuff circle
circle = Circle(
    (nerve_cx_px, nerve_cy_px),
    subject_radius_px,
    fill=False,
    edgecolor="cyan",
    linewidth=2
)
plt.gca().add_patch(circle)

plt.title(f"Axon centers + stim cuff circle for {subject_to_plot}")
plt.axis("off")
plt.show()

# ---------------------------------------------------------
# Plot 2: add a line showing shortest distance for one axon
# ---------------------------------------------------------
example_row = sub_hist.iloc[0]
ax_x_px = example_row["CentroidX_px"]
ax_y_px = example_row["CentroidY_px"]

plt.figure(figsize=(8, 8))
plt.imshow(img, cmap="gray")
plt.scatter(sub_hist["CentroidX_px"], sub_hist["CentroidY_px"], s=5, c="red")
circle = Circle(
    (nerve_cx_px, nerve_cy_px),
    subject_radius_px,
    fill=False,
    edgecolor="cyan",
    linewidth=2
)
plt.gca().add_patch(circle)

# Compute direction vector from nerve center → axon
d_px = np.sqrt((ax_x_px - nerve_cx_px)**2 + (ax_y_px - nerve_cy_px)**2)
ux = (ax_x_px - nerve_cx_px) / d_px
uy = (ax_y_px - nerve_cy_px) / d_px

# Boundary point on cuff circle
boundary_x = nerve_cx_px + ux * subject_radius_px
boundary_y = nerve_cy_px + uy * subject_radius_px

# Draw distance line
plt.plot([ax_x_px, boundary_x], [ax_y_px, boundary_y], "yellow")

plt.title(f"Example distance measurement for {subject_to_plot}")
plt.axis("off")
plt.show()
