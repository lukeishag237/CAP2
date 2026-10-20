import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import yaml

# ============================================================
# Load config and CNAP dataframe
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = local_config["local_output_path"]
RAW_DATA_PATH = local_config["raw_data_path"]

# STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")
# cnap_df = pd.read_pickle(STORED_DF_PATH)
STORED_DF_PATH  = os.path.join(LOCAL_OUTPUT_PATH, "stored_df_edited.csv")
cnap_df = pd.read_csv(STORED_DF_PATH)

# ============================================================
# Load histology data (NEW pipeline output)
# ============================================================
HIST_CSV = os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv")
hist_df = pd.read_csv(HIST_CSV)

# New column names from your updated pipeline:
#   - MinFeret
#   - distance_to_stim_cuff_um
#   - NerveType
#   - Name

nerve_map = {"sciatic": "SN", "tibial": "tibial", "vagus": "cVN"}
hist_df["NerveCanonical"] = hist_df["NerveType"].map(nerve_map)

# Expdate mapping (same as before)
mapping = {
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
hist_df["expdate"] = hist_df["Name"].map(mapping)

# Drop rows without mapping or canonical nerve
hist_df = hist_df.dropna(subset=["expdate", "NerveCanonical"])

# Ensure expdate is string
hist_df["expdate"] = hist_df["expdate"].astype(str)
cnap_df["expdate"] = cnap_df["expdate"].astype(str)

print("Histology columns:", list(hist_df.columns))

# ============================================================
# Aggregate histology → nerve-level summaries (NEW: MinFeret)
# ============================================================
hist_summary = (
    hist_df.groupby(["expdate", "NerveCanonical"], as_index=False)
           .agg(
               hist_mean_diam=("MinFeret", "mean"),
               hist_median_diam=("MinFeret", "median"),
               hist_total_diam=("MinFeret", "sum"),
               hist_axon_count=("MinFeret", "count")
           )
)

# ============================================================
# Basic CNAP sanity checks
# ============================================================
print("\nCNAP peak_time summary:")
print(cnap_df["peak_time"].describe())

print("\nCNAP stim_rec_dist summary:")
print(cnap_df["stim_rec_dist"].describe())

# Keep only rows with valid latency and distance
cnap_valid = cnap_df.dropna(subset=["peak_time", "stim_rec_dist"]).copy()
cnap_valid = cnap_valid[(cnap_valid["peak_time"] > 0) & (cnap_valid["stim_rec_dist"] > 0)]

# ============================================================
# Compute conduction velocity (m/s)
# ============================================================
# Fix pathological outliers (same logic as before)
for _ in range(3):
    idx = cnap_valid["peak_time"].idxmin()
    cnap_valid.loc[idx, "peak_time"] = 0.6

cnap_valid["cv_m_per_s"] = (
    cnap_valid["stim_rec_dist"] / cnap_valid["peak_time"]
)

print("\nConduction velocity summary (all CAPs):")
print(cnap_valid["cv_m_per_s"].describe())

# ============================================================
# Focus on myelinated fibers and nerves of interest
# ============================================================
NERVES_OF_INTEREST = ["SN", "tibial", "cVN"]

cnap_myel = cnap_valid[
    (cnap_valid["target_fiber_types"] == "myel") &
    (cnap_valid["rec_nerve"].isin(NERVES_OF_INTEREST))
].copy()

print("\nCounts by nerve (myelinated only):")
print(cnap_myel["rec_nerve"].value_counts())

# ============================================================
# Merge CNAP with histology summaries
# ============================================================
merged = pd.merge(
    cnap_myel,
    hist_summary,
    left_on=["expdate", "rec_nerve"],
    right_on=["expdate", "NerveCanonical"],
    how="left"
)

print("\nMerged CNAP + histology shape:", merged.shape)

# ============================================================
# Compare CV across nerves
# ============================================================
summary_by_nerve = (
    merged.groupby("rec_nerve")["cv_m_per_s"]
          .agg(["count", "mean", "std", "median"])
          .reset_index()
)

print("\nConduction velocity by nerve (myelinated only):")
print(summary_by_nerve)

# ============================================================
# CV distribution plot (unchanged)
# ============================================================
plt.figure(figsize=(8, 6))
sns.boxplot(
    data=merged,
    x="rec_nerve",
    y="cv_m_per_s",
    order=["cVN", "SN", "tibial"]
)
sns.stripplot(
    data=merged,
    x="rec_nerve",
    y="cv_m_per_s",
    order=["cVN", "SN", "tibial"],
    color="black",
    alpha=0.4
)
plt.xlabel("Nerve")
plt.ylabel("Conduction velocity (m/s)")
plt.title("Conduction velocity by nerve (myelinated CAPs)")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# ============================================================
# Diameter-normalized CV analysis
# ============================================================

print("\n=== Diameter-Normalized CV Analysis ===\n")

# 1. Predict CV per axon using MinFeret
k_factor = 5.7
hist_df["predicted_cv_m_per_s"] = hist_df["MinFeret"] * k_factor

pred_summary = (
    hist_df.groupby(["expdate", "NerveCanonical"], as_index=False)
           .agg(predicted_cv_mean=("predicted_cv_m_per_s", "mean"))
)

# 2. Merge predicted CV with measured CV
merged_pred = pd.merge(
    merged,
    pred_summary,
    left_on=["expdate", "rec_nerve"],
    right_on=["expdate", "NerveCanonical"],
    how="left"
)

print(merged_pred[["expdate", "rec_nerve", "predicted_cv_mean"]].head(20))

# 3. Plot predicted vs measured CV (unchanged)
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=merged_pred,
    x="predicted_cv_mean",
    y="cv_m_per_s",
    hue="rec_nerve",
    alpha=0.7,
    s=80
)

max_val = max(merged_pred["predicted_cv_mean"].max(),
              merged_pred["cv_m_per_s"].max())
plt.plot([0, max_val], [0, max_val], 'k--')

plt.xlabel("Predicted CV from diameter (m/s)")
plt.ylabel("Measured CAP CV (m/s)")
plt.title("Predicted vs Measured Conduction Velocity")
plt.legend(title="Nerve")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# 4. Residuals
merged_pred["cv_residual"] = (
    merged_pred["cv_m_per_s"] - merged_pred["predicted_cv_mean"]
)

# 5. Residual plot (unchanged)
plt.figure(figsize=(8, 6))
sns.boxplot(
    data=merged_pred,
    x="rec_nerve",
    y="cv_residual",
    order=["cVN", "SN", "tibial"]
)
sns.stripplot(
    data=merged_pred,
    x="rec_nerve",
    y="cv_residual",
    order=["cVN", "SN", "tibial"],
    color="black",
    alpha=0.4
)

plt.axhline(0, color='k', linestyle='--')
plt.xlabel("Nerve")
plt.ylabel("Residual CV (Measured – Predicted)")
plt.title("Diameter-Normalized Conduction Velocity Differences")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.show()

# 6. Mean diameter vs measured CV (unchanged)
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=merged_pred,
    x="hist_mean_diam",
    y="cv_m_per_s",
    hue="rec_nerve",
    alpha=0.7,
    s=80
)

plt.xlabel("Mean axon diameter (µm)")
plt.ylabel("Measured CAP CV (m/s)")
plt.title("Mean Axon Diameter vs Measured CV")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("\n=== Diameter-normalized CV analysis complete ===\n")
