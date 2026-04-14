import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import yaml
import os

# ============================================================
# Load configuration and CNAP data
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = local_config["local_output_path"]
OUTPUT_DIR = LOCAL_OUTPUT_PATH + "/OutputFigs"
# STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")
# cnap_df = pd.read_pickle(STORED_DF_PATH)
STORED_DF_PATH  = os.path.join(LOCAL_OUTPUT_PATH, "stored_df_edited.csv")
cnap_df = pd.read_csv(STORED_DF_PATH)


# ============================================================
# Load histology data (axon-level)
# ============================================================
HIST_CSV = os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv")
hist_df = pd.read_csv(HIST_CSV)

nerve_map = {"sciatic": "SN", "tibial": "tibial", "vagus": "cVN"}

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
hist_df["NerveCanonical"] = hist_df["NerveType"].map(nerve_map)
hist_df = hist_df.dropna(subset=["expdate", "NerveCanonical"])

hist_df["expdate"] = hist_df["expdate"].astype(str)
cnap_df["expdate"] = cnap_df["expdate"].astype(str)


# ============================================================
# Aggregate axon-level histology → nerve-level summaries
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
# Load fascicle-level data
# ============================================================
RAW_DATA_PATH = local_config.get("raw_data_path")
BASE_SUMMARY_DIR = os.path.join(RAW_DATA_PATH, "..", "60x", "derivatives", "summary")

FASCICLE_CSV = os.path.join(BASE_SUMMARY_DIR, "fascicle_data.csv")
fas_df = pd.read_csv(FASCICLE_CSV)

fas_df["expdate"] = fas_df["subject"].map(mapping)
fas_df["expdate"] = fas_df["expdate"].astype(str)


nerve_lookup = hist_df[["Name", "NerveType"]].drop_duplicates()
fas_df = fas_df.merge(nerve_lookup, left_on="subject", right_on="Name", how="left")
fas_df["NerveCanonical"] = fas_df["NerveType"].map(nerve_map)

# ============================================================
# Extract summary rows (nerve-level fascicle metrics)
# ============================================================
fas_summary_rows = fas_df[fas_df["mask_type"] == "summary"].copy()

fas_summary = fas_summary_rows[[
    "expdate", "NerveCanonical",
    "num_fascicles",
    "total_area_um2",
    "axon_area_um2",
    "area_ratio",
    "outer_area_um2",
    "outer_perimeter_um",
    "outer_thickness_um",
    "whole_area_um2"
]]

# Compute area proportion
fas_summary["area_proportion"] = (
    fas_summary["total_area_um2"] / fas_summary["outer_area_um2"]
)

# ============================================================
# Compute mean distance_to_electrode per subject (axon → nerve)
# ============================================================
dist_summary = (
    hist_df.groupby(["Name", "expdate"], as_index=False)
           .agg(mean_distance_to_electrode=("distance_to_electrode", "mean"))
)

# ============================================================
# Merge CNAP + histology + fascicle
# ============================================================
NERVES = ["cVN", "SN", "tibial"]

df = cnap_df[
    (cnap_df["rec_nerve"].isin(NERVES)) &
    (cnap_df["ch_num_neural"] == 1) &
    (cnap_df["stim_rec_dist"] > 0)
].copy()

df = pd.merge(
    df,
    hist_summary,
    left_on=["expdate", "rec_nerve"],
    right_on=["expdate", "NerveCanonical"],
    how="inner"
)

df = pd.merge(
    df,
    fas_summary,
    on=["expdate", "NerveCanonical"],
    how="inner"
)

df = df.merge(
    dist_summary[["expdate", "mean_distance_to_electrode"]],
    on="expdate",
    how="left"
).rename(columns={"mean_distance_to_electrode": "distance_to_electrode"})

# ============================================================
# Prepare variables
# ============================================================
df = df[df["peak_amplitude"] > 0]
df["fiber"] = df["target_fiber_types"].astype("category")
df["nerve"] = df["rec_nerve"].astype("category")
df["subject"] = df["expdate"].astype("category")

# ============================================================
# Classify cVN myelinated vs unmyelinated
# ============================================================
def classify_cvn(row):
    if row["nerve"] != "cVN":
        return row["nerve"]
    return "cVN_myel" if row["fiber"] == "myel" else "cVN_unmyel"

df["nerve_fiber_group"] = df.apply(classify_cvn, axis=1)

nerve_groups = {
    "tibial": "tibial",
    "sciatic": "SN",
    "cVN_myel": "cVN_myel",
    #"cVN_unmyel": "cVN_unmyel"
}

# ============================================================
# Compute average fascicle size
# ============================================================
df["avg_fascicle_size"] = df["total_area_um2"] / df["num_fascicles"]

# ============================================================
# Stage 2A: Multi-group direct PEAK AMPLITUDE models
# ============================================================
PREDICTORS = {
    "stim_rec_dist": "Conduction Distance (mm)",
    "hist_axon_count": "Axon Count",
    "total_area_um2": "Inner Area (µm²)",
    "area_ratio": "Area Ratio (f/i)",
    "avg_fascicle_size": "Average Fascicle Size (µm²)",
    "outer_thickness_um": "Outer Thickness (µm)",
    "area_proportion": "Area Proportion (inner/outer)",
    "distance_to_electrode": "Electrode Fit (µm)",
}

colors = {
    "tibial": "tab:blue",
    "sciatic": "tab:green",
    "cVN_myel": "tab:red"
}

print("\n================ DIRECT PEAK AMPLITUDE MODELS ================\n")

for pred, pred_label in PREDICTORS.items():
    print(f"\n========== PREDICTOR: {pred_label} ==========\n")

    plt.figure(figsize=(8, 7))

    for label, code in nerve_groups.items():
        df_sub = df[df["nerve_fiber_group"] == code].copy()
        if df_sub.empty:
            continue

        df_valid = df_sub.dropna(subset=["peak_amplitude", pred])
        if df_valid.empty:
            continue

        model_uni = smf.ols(f"peak_amplitude ~ {pred}", data=df_valid).fit()
        r2 = model_uni.rsquared
        pval = model_uni.pvalues[pred]

        plt.scatter(
            df_valid[pred],
            df_valid["peak_amplitude"],
            alpha=0.6,
            color=colors[label]
        )

        x_vals = np.linspace(df_valid[pred].min(), df_valid[pred].max(), 200)
        y_hat = model_uni.params["Intercept"] + model_uni.params[pred] * x_vals

        plt.plot(
            x_vals,
            y_hat,
            color=colors[label],
            linewidth=2,
            label=f"{label} (R²={r2:.2f}, p={pval:.3f})"
        )

    plt.xlabel(pred_label)
    plt.ylabel("Peak Amplitude")
    plt.title(f"Peak Amplitude vs {pred_label} — All Nerves")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # SAVE FIGURE
    fig_path = os.path.join(OUTPUT_DIR, f"PeakAmplitude_vs_{pred}.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()


# ============================================================
# Stage 2B: Multi-group direct PEAK TIME models
# ============================================================
print("\n================ DIRECT PEAK TIME MODELS ================\n")

for pred, pred_label in PREDICTORS.items():
    print(f"\n========== PREDICTOR: {pred_label} ==========\n")

    plt.figure(figsize=(8, 7))

    for label, code in nerve_groups.items():
        df_sub = df[df["nerve_fiber_group"] == code].copy()
        if df_sub.empty:
            continue

        df_valid = df_sub.dropna(subset=["peak_time", pred])
        if df_valid.empty:
            continue

        model_uni = smf.ols(f"peak_time ~ {pred}", data=df_valid).fit()
        r2 = model_uni.rsquared
        pval = model_uni.pvalues[pred]

        plt.scatter(
            df_valid[pred],
            df_valid["peak_time"],
            alpha=0.6,
            color=colors[label]
        )

        x_vals = np.linspace(df_valid[pred].min(), df_valid[pred].max(), 200)
        y_hat = model_uni.params["Intercept"] + model_uni.params[pred] * x_vals

        plt.plot(
            x_vals,
            y_hat,
            color=colors[label],
            linewidth=2,
            label=f"{label} (R²={r2:.2f}, p={pval:.3f})"
        )

    plt.xlabel(pred_label)
    plt.ylabel("Peak Time (ms)")
    plt.title(f"Peak Time vs {pred_label} — All Nerves")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # SAVE FIGURE
    fig_path = os.path.join(OUTPUT_DIR, f"PeakTime_vs_{pred}.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
