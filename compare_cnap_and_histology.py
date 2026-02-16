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
STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")
cnap_df = pd.read_pickle(STORED_DF_PATH)

# ============================================================
# Load histology data (axon-level)
# ============================================================
RAW_DATA_PATH = local_config.get("raw_data_path")
BASE_SUMMARY_DIR = os.path.join(RAW_DATA_PATH, "..", "60x", "derivatives", "summary")

HIST_CSV = os.path.join(BASE_SUMMARY_DIR, "all_samples_fiber_diameters.csv")
hist_df = pd.read_csv(HIST_CSV)

nerve_map = {"sciatic": "SN", "tibial": "tibial", "vagus": "cVN"}
hist_df["NerveCanonical"] = hist_df["NerveType"].map(nerve_map)

mapping = {
    "39-1": "20220531", "39-2": "20220531",
    "40-1": "20220622",
    "41-1": "20220720", "41-2": "20220720",
    "43-1": "20221103",
    "44-1": "20221108", "44-2": "20221108",
    "45-1": "20221122",
    "46-1": "20221201", "46-2": "20221202",
    "47-1": "20230131",
    "48-1": "20230207", "48-2": "20230207"
}
hist_df["expdate"] = hist_df["Name"].map(mapping)
hist_df = hist_df.dropna(subset=["expdate", "NerveCanonical"])

# ============================================================
# Aggregate axon-level histology → nerve-level summaries
# ============================================================
hist_summary = (
    hist_df.groupby(["expdate", "NerveCanonical"], as_index=False)
           .agg(
               hist_mean_diam=("SAE_diameter", "mean"),
               hist_median_diam=("SAE_diameter", "median"),
               hist_total_diam=("SAE_diameter", "sum"),
               hist_axon_count=("SAE_diameter", "count")
           )
)

# ============================================================
# Load fascicle-level data (from fasicle_data.csv)
# ============================================================
FASCICLE_CSV = os.path.join(BASE_SUMMARY_DIR, "fasicle_data.csv")
fas_df = pd.read_csv(FASCICLE_CSV)

# 1. Map subject → expdate using the same mapping as histology
fas_df["expdate"] = fas_df["subject"].map(mapping)

# 2. Merge in nerve type from histology (Name == subject)
nerve_lookup = hist_df[["Name", "NerveType"]].drop_duplicates()
fas_df = fas_df.merge(nerve_lookup, left_on="subject", right_on="Name", how="left")

# 3. Map NerveType → NerveCanonical
if "NerveType" in fas_df.columns:
    fas_df["NerveCanonical"] = fas_df["NerveType"].map(nerve_map)

# Split inner (i) and outer (o) metrics and aggregate to nerve-level
fas_i = fas_df[fas_df["mask_type"] == "i"].copy()
fas_o = fas_df[fas_df["mask_type"] == "o"].copy()

fas_i_summary = (
    fas_i.groupby(["expdate", "NerveCanonical"], as_index=False)
         .agg(
             i_eccentricity=("eccentricity", "mean"),
             i_area=("area", "mean"),
             i_orientation=("orientation", "mean")
         )
)

fas_o_summary = (
    fas_o.groupby(["expdate", "NerveCanonical"], as_index=False)
         .agg(
             o_thickness=("thickness_estimate", "mean"),
             o_area=("area", "mean")
         )
)

fas_summary = pd.merge(
    fas_i_summary,
    fas_o_summary,
    on=["expdate", "NerveCanonical"],
    how="outer"
)

# ============================================================
# Merge histology + fascicle summaries
# ============================================================
hist_fas_summary = pd.merge(
    hist_summary,
    fas_summary,
    on=["expdate", "NerveCanonical"],
    how="inner"
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
    hist_fas_summary,
    left_on=["expdate", "rec_nerve"],
    right_on=["expdate", "NerveCanonical"],
    how="inner"
)

# ============================================================
# Prepare variables (ONLY peak_amplitude)
# ============================================================
df = df[df["peak_amplitude"] > 0]
df["log_peak_amp"] = np.log(df["peak_amplitude"])

df["log_dist"] = np.log(df["stim_rec_dist"])
df["fiber"] = df["target_fiber_types"].astype("category")
df["nerve"] = df["rec_nerve"].astype("category")
df["subject"] = df["expdate"].astype("category")

# ============================================================
# Split cVN into myelinated vs unmyelinated
# ============================================================
def classify_cvn(row):
    if row["nerve"] != "cVN":
        return row["nerve"]  # SN or tibial
    if row["fiber"] == "myel":
        return "cVN_myel"
    elif row["fiber"] == "unmyel":
        return "cVN_unmyel"
    return "cVN"

df["nerve_fiber_group"] = df.apply(classify_cvn, axis=1)
df["nerve_fiber_group"] = df["nerve_fiber_group"].astype("category")

# ============================================================
# Helper: distance-only OLS model log(peak_amp) ~ log(distance)
# ============================================================
def fit_distance_only(df_sub):
    model = smf.ols("log_peak_amp ~ log_dist", data=df_sub).fit()
    df_sub["log_peak_amp_dist_resid"] = model.resid
    return model, df_sub

# ============================================================
# Nerve groups
# ============================================================
nerve_groups = {
    "tibial": "tibial",
    "sciatic": "SN",
    "cVN_myel": "cVN_myel",
    "cVN_unmyel": "cVN_unmyel"
}

# Predictors for second-stage models
SECOND_STAGE_PREDICTORS = [
    "hist_axon_count",
    "hist_mean_diam",
    "i_eccentricity",
    "i_area",
    "o_thickness",
    "i_orientation",
]

# ============================================================
# Distance-only model per nerve-fiber group
# + Option A only: univariate models + plots (24 total)
# ============================================================
print("\n================ METRIC: PEAK_AMPLITUDE =================\n")

for label, code in nerve_groups.items():
    print(f"\n---- GROUP: {label.upper()} ----")

    df_sub = df[df["nerve_fiber_group"] == code].copy()
    if df_sub.empty:
        print("No data for this group.")
        continue

    # Distance-only OLS model
    result_ols, df_sub = fit_distance_only(df_sub)
    print("\nDistance-only OLS model summary (log_peak_amp ~ log_dist):")
    print(result_ols.summary())

    resid_col = "log_peak_amp_dist_resid"

    # ----------------------------------------------------
    # Option A: Univariate OLS models (resid ~ predictor)
    # + individual plots
    # ----------------------------------------------------
    print("\nOption A: Univariate second-stage models (resid ~ predictor)\n")
    for pred in SECOND_STAGE_PREDICTORS:
        if pred not in df_sub.columns:
            print(f"  Skipping {pred}: not in dataframe.")
            continue

        df_valid = df_sub.dropna(subset=[resid_col, pred])
        if df_valid.empty:
            print(f"  Skipping {pred}: no valid data.")
            continue

        formula = f"{resid_col} ~ {pred}"
        model_uni = smf.ols(formula, data=df_valid).fit()
        print(f"\n  Predictor: {pred}")
        print(model_uni.summary())

        # Plot residuals vs this predictor
        plt.figure(figsize=(7, 6))
        plt.scatter(
            df_valid[pred],
            df_valid[resid_col],
            alpha=0.6
        )
        # Add fitted line
        x_vals = np.linspace(df_valid[pred].min(), df_valid[pred].max(), 100)
        y_hat = model_uni.params["Intercept"] + model_uni.params[pred] * x_vals
        plt.plot(x_vals, y_hat, color="red", linewidth=2)

        plt.axhline(0, color="black", linestyle="--")
        plt.xlabel(pred)
        plt.ylabel("Peak amplitude residual (after distance)")
        plt.title(f"PEAK_AMPLITUDE — {label} — resid vs {pred}")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
