import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io

# ---------------------------------------------------------
# Load config and paths
# ---------------------------------------------------------
CONFIG_FILE = "local_config.yml"
with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

BASE_DIR = os.path.join(
    local_config["raw_data_path"],
    "../60x/derivatives/segmentations/manually_curated"
)

CSV_PATH = os.path.join(
    local_config["raw_data_path"],
    "..", "60x", "derivatives", "summary", "fascicle_data.csv"
)

df = pd.read_csv(CSV_PATH)

# Keep only fascicle rows (mask_type == "i")
fas_df = df[df["mask_type"] == "i"].copy()

# ---------------------------------------------------------
# Helper: draw ellipse from regionprops metrics
# ---------------------------------------------------------
def plot_ellipse(ax, centroid_row, centroid_col, major, minor, orientation, color="red"):
    """
    Draw an ellipse using regionprops-style parameters.
    major/minor are full lengths, so radius = length/2.
    orientation is in radians, CCW from horizontal axis.
    """
    t = np.linspace(0, 2*np.pi, 200)
    a = major / 2
    b = minor / 2

    # Parametric ellipse before rotation
    x = a * np.cos(t)
    y = b * np.sin(t)

    # Apply rotation
    theta = -orientation + np.pi/2
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    xy = np.vstack((x, y))
    xy_rot = R @ xy

    # Translate to centroid
    x_final = xy_rot[0, :] + centroid_col
    y_final = xy_rot[1, :] + centroid_row

    ax.plot(x_final, y_final, color=color, linewidth=2)


# ---------------------------------------------------------
# MAIN LOOP: plot all fascicle ellipses over i.tif
# ---------------------------------------------------------
for subject, sub_df in fas_df.groupby("subject"):

    print(f"\nProcessing subject: {subject}")

    subject_dir = os.path.join(BASE_DIR, subject)
    if not os.path.exists(subject_dir):
        print(f"  WARNING: directory not found: {subject_dir}")
        continue

    # -----------------------------------------------------
    # Load the canonical background mask: i.tif
    # -----------------------------------------------------
    background_path = os.path.join(subject_dir, "i.tif")
    if not os.path.exists(background_path):
        print(f"  ERROR: i.tif not found in {subject_dir}")
        continue

    background_mask = io.imread(background_path)

    # -----------------------------------------------------
    # Plot background once, overlay all ellipses
    # -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(background_mask, cmap="gray")
    ax.set_title(f"{subject} — All fascicle ellipses")

    # Overlay each fascicle ellipse
    for _, row in sub_df.iterrows():
        plot_ellipse(
            ax,
            centroid_row=row["centroid_row"],
            centroid_col=row["centroid_col"],
            major=row["major_axis_length"],
            minor=row["minor_axis_length"],
            orientation=row["orientation"],
            color="red"
        )

    ax.set_axis_off()
    plt.tight_layout()
    plt.show()
