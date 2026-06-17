"""
fit_distance_models.py

Two-stage modelling of CNAP features vs. histology/morphology predictors.

Stage 1 — Distance-only baseline:
    peak_amplitude ~ stim_rec_dist
    peak_time      ~ stim_rec_dist

Stage 2 — Distance-adjusted models for each additional predictor:
    peak_amplitude ~ stim_rec_dist + <predictor>
    peak_time      ~ stim_rec_dist + <predictor>

All figures are shown (not saved).
"""

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
RAW_DATA_PATH     = local_config.get("raw_data_path")

STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")
cnap_df = pd.read_pickle(STORED_DF_PATH)

# ============================================================
# Load histology data (axon-level)
# ============================================================
HIST_CSV = os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv")
hist_df  = pd.read_csv(HIST_CSV)

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

hist_df["expdate"]        = hist_df["Name"].map(mapping)
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
               hist_mean_diam  = ("MinFeret", "mean"),
               hist_median_diam= ("MinFeret", "median"),
               hist_total_diam = ("MinFeret", "sum"),
               hist_axon_count = ("MinFeret", "count"),
           )
)

# ============================================================
# Load fascicle-level data
# ============================================================
BASE_SUMMARY_DIR = os.path.join(RAW_DATA_PATH, "..", "60x", "derivatives", "summary")
fas_df = pd.read_csv(os.path.join(BASE_SUMMARY_DIR, "fascicle_data.csv"))

fas_df["expdate"] = fas_df["subject"].map(mapping)
fas_df["expdate"] = fas_df["expdate"].astype(str)

nerve_lookup = hist_df[["Name", "NerveType"]].drop_duplicates()
fas_df = fas_df.merge(nerve_lookup, left_on="subject", right_on="Name", how="left")
fas_df["NerveCanonical"] = fas_df["NerveType"].map(nerve_map)

fas_summary = (
    fas_df[fas_df["mask_type"] == "summary"]
    [[
        "expdate", "NerveCanonical",
        "num_fascicles", "total_area_um2", "axon_area_um2",
        "area_ratio", "outer_area_um2", "outer_perimeter_um",
        "outer_thickness_um", "whole_area_um2",
    ]]
    .copy()
)
fas_summary["area_proportion"] = (
    fas_summary["total_area_um2"] / fas_summary["outer_area_um2"]
)

# ============================================================
# Mean distance-to-electrode per subject
# ============================================================
dist_summary = (
    hist_df.groupby(["Name", "expdate"], as_index=False)
           .agg(mean_distance_to_electrode=("distance_to_stim_cuff_um", "mean"))
)

# ============================================================
# Merge CNAP + histology + fascicle
# ============================================================
NERVES = ["cVN", "SN", "tibial"]

df = cnap_df[
    (cnap_df["rec_nerve"].isin(NERVES)) &
    (cnap_df["stim_rec_dist"] > 0)
].copy()

# Join dataset_electrodetype from Excel metadata so we know how many
# contacts each recording has (2 or 3)
_meta = pd.read_excel(os.path.join(RAW_DATA_PATH, "AnalysisScript_CAPs_Batch_parameters.xlsx"))
_meta["expdate"] = _meta["file_prefix_outputs_values"].str[:8]
_elec = (
    _meta.groupby(["expdate", "rec_nerve"], as_index=False)
         .agg(dataset_electrodetype=("dataset_electrodetype", "first"))
)
df = df.merge(_elec, on=["expdate", "rec_nerve"], how="left")

# Apply per-channel distance offset based on electrode geometry:
#   2-contact cuff: ch1 = -0.75 mm, ch2 = +0.75 mm
#   3-contact cuff: ch1 = -1.5 mm,  ch2 = 0,       ch3 = +1.5 mm
def channel_offset(row):
    n = row["dataset_electrodetype"]
    ch = row["ch_num_neural"]
    if n == 2:
        offsets = {1: -0.75, 2: +0.75}
    elif n == 3:
        offsets = {1: -1.5, 2: 0.0, 3: +1.5}
    else:
        offsets = {}
    return offsets.get(ch, 0.0)

df["stim_rec_dist"] = df["stim_rec_dist"] + df.apply(channel_offset, axis=1)

df = pd.merge(df, hist_summary,
              left_on=["expdate", "rec_nerve"],
              right_on=["expdate", "NerveCanonical"],
              how="inner")

df = pd.merge(df, fas_summary,
              on=["expdate", "NerveCanonical"],
              how="inner")

df = df.merge(
    dist_summary[["expdate", "mean_distance_to_electrode"]],
    on="expdate", how="left"
).rename(columns={"mean_distance_to_electrode": "distance_to_stim_cuff_um"})

# ============================================================
# Load impedance data and join per channel + nerve
# ============================================================
imp_raw = pd.read_csv("all_impedences.csv")
imp_raw["expdate"] = imp_raw["expdate"].astype(str)

# Build a lookup: Histology ID → NerveCanonical using the same mapping + nerve_map
_hist_nerve = (
    hist_df[["Name", "NerveCanonical", "expdate"]]
    .drop_duplicates()
    .rename(columns={"Name": "Histology ID"})
)
imp_raw = imp_raw.merge(_hist_nerve, on=["Histology ID", "expdate"], how="left")

# Melt ch_1_vs_ref / ch_2_vs_ref / ch_3_vs_ref → long format keyed on
# expdate + NerveCanonical + ch_num_neural
imp_long = pd.melt(
    imp_raw,
    id_vars=["expdate", "NerveCanonical"],
    value_vars=["ch_1_vs_ref", "ch_2_vs_ref", "ch_3_vs_ref"],
    var_name="ch_col",
    value_name="impedance_kOhm"
)
imp_long["ch_num_neural"] = imp_long["ch_col"].map(
    {"ch_1_vs_ref": 1, "ch_2_vs_ref": 2, "ch_3_vs_ref": 3}
)
imp_long = imp_long.drop(columns="ch_col").dropna(subset=["NerveCanonical"])

df = df.merge(
    imp_long[["expdate", "NerveCanonical", "ch_num_neural", "impedance_kOhm"]],
    left_on=["expdate", "NerveCanonical", "ch_num_neural"],
    right_on=["expdate", "NerveCanonical", "ch_num_neural"],
    how="left"
)

# ============================================================
# Prepare variables
# ============================================================
df = df[df["peak_amplitude"] > 0].copy()
df["fiber"]   = df["target_fiber_types"].astype("category")
df["nerve"]   = df["rec_nerve"].astype("category")
df["subject"] = df["expdate"].astype("category")
df["avg_fascicle_size"] = df["total_area_um2"] / df["num_fascicles"]

def classify_cvn(row):
    if row["nerve"] != "cVN":
        return row["nerve"]
    return "cVN_myel" if row["fiber"] == "myel" else "cVN_unmyel"

df["nerve_fiber_group"] = df.apply(classify_cvn, axis=1)

nerve_groups = {
    "tibial":   "tibial",
    "SN":       "SN",
    "cVN_myel": "cVN_myel",
}

colors = {
    "tibial":   "tab:blue",
    "SN":       "tab:green",
    "cVN_myel": "tab:red",
}

# Non-distance predictors (used in Stage 2 & 4 residual models)
# impedance_kOhm is handled separately in Stage 3
EXTRA_PREDICTORS = {
    "hist_axon_count":          "Axon Count",
    "total_area_um2":           "Inner Area (µm²)",
    "area_ratio":               "Area Ratio (f/i)",
    "avg_fascicle_size":        "Average Fascicle Size (µm²)",
    "outer_thickness_um":       "Outer Thickness (µm)",
    "area_proportion":          "Area Proportion (inner/outer)",
    "distance_to_stim_cuff_um": "Electrode Fit (µm)",
}

OUTCOMES = {
    "peak_amplitude": "Peak Amplitude (µV)",
    "peak_time":      "Peak Time (ms)",
}

# ============================================================
# Stage 1 — Distance-only baseline
# ============================================================
print("\n================ STAGE 1: DISTANCE BASELINE ================\n")

for outcome, outcome_label in OUTCOMES.items():
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, code in nerve_groups.items():
        df_sub   = df[df["nerve_fiber_group"] == code].copy()
        df_valid = df_sub.dropna(subset=[outcome, "stim_rec_dist"])
        if df_valid.empty:
            continue
        model = smf.ols(f"{outcome} ~ stim_rec_dist", data=df_valid).fit()
        print(f"  [{label}] {outcome} ~ stim_rec_dist  "
              f"R²={model.rsquared:.3f}  p={model.pvalues['stim_rec_dist']:.4f}")
        ax.scatter(df_valid["stim_rec_dist"], df_valid[outcome],
                   alpha=0.6, color=colors[label], zorder=3)
        x_vals = np.linspace(df_valid["stim_rec_dist"].min(),
                             df_valid["stim_rec_dist"].max(), 200)
        y_hat = (model.params["Intercept"]
                 + model.params["stim_rec_dist"] * x_vals)
        ax.plot(x_vals, y_hat, color=colors[label], linewidth=2,
                label=f"{label} (R²={model.rsquared:.2f}, "
                      f"p={model.pvalues['stim_rec_dist']:.3f})")
    ax.set_xlabel("Conduction Distance (mm)", fontsize=14)
    ax.set_ylabel(outcome_label, fontsize=14)
    ax.set_title(f"{outcome_label} ~ Conduction Distance", fontsize=15)
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

# ============================================================
# Stage 1b — Impedance-only baseline
# ============================================================
print("\n================ STAGE 1b: IMPEDANCE BASELINE ================\n")

for outcome, outcome_label in OUTCOMES.items():
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, code in nerve_groups.items():
        df_sub   = df[df["nerve_fiber_group"] == code].copy()
        df_valid = df_sub.dropna(subset=[outcome, "impedance_kOhm"])
        if df_valid.empty:
            continue
        model = smf.ols(f"{outcome} ~ impedance_kOhm", data=df_valid).fit()
        print(f"  [{label}] {outcome} ~ impedance_kOhm  "
              f"R²={model.rsquared:.3f}  p={model.pvalues['impedance_kOhm']:.4f}")
        ax.scatter(df_valid["impedance_kOhm"], df_valid[outcome],
                   alpha=0.6, color=colors[label], zorder=3)
        x_vals = np.linspace(df_valid["impedance_kOhm"].min(),
                             df_valid["impedance_kOhm"].max(), 200)
        y_hat = (model.params["Intercept"]
                 + model.params["impedance_kOhm"] * x_vals)
        ax.plot(x_vals, y_hat, color=colors[label], linewidth=2,
                label=f"{label} (R²={model.rsquared:.2f}, "
                      f"p={model.pvalues['impedance_kOhm']:.3f})")
    ax.set_xlabel("Impedance (kΩ)", fontsize=14)
    ax.set_ylabel(outcome_label, fontsize=14)
    ax.set_title(f"{outcome_label} ~ Impedance", fontsize=15)
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

# ============================================================
# Stage 2 — Dist+impedance residuals vs each predictor
# Two figures per predictor (peak_amplitude then peak_time).
# Each figure is 1×2: left = raw scatter (no model), right = residuals
# ============================================================
print("\n================ STAGE 2: DIST+IMP RESIDUALS ================\n")

for pred, pred_label in EXTRA_PREDICTORS.items():
    print(f"\n  Predictor: {pred}")

    for outcome, outcome_label in OUTCOMES.items():
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        ax_raw   = axes[0]
        ax_resid = axes[1]

        for label, code in nerve_groups.items():
            df_sub   = df[df["nerve_fiber_group"] == code].copy()
            df_valid = df_sub.dropna(
                subset=[outcome, "stim_rec_dist", "impedance_kOhm", pred])
            if len(df_valid) < 3:
                continue

            # Dist + impedance model → residuals
            base_model  = smf.ols(
                f"{outcome} ~ stim_rec_dist + impedance_kOhm", data=df_valid).fit()
            df_valid    = df_valid.copy()
            df_valid["resid"] = base_model.resid
            resid_model = smf.ols("resid ~ " + pred, data=df_valid).fit()
            r2_resid    = resid_model.rsquared
            p_resid     = resid_model.pvalues.get(pred, np.nan)

            print(f"    [{outcome}] {label}:  "
                  f"R²(resid)={r2_resid:.3f}  p={p_resid:.4f}")

            x_p = np.linspace(df_valid[pred].min(), df_valid[pred].max(), 200)

            # Left: raw scatter only, no model line
            ax_raw.scatter(df_valid[pred], df_valid[outcome],
                           alpha=0.6, color=colors[label], zorder=3,
                           label=label)

            # Right: dist+imp residuals + trendline
            ax_resid.scatter(df_valid[pred], df_valid["resid"],
                             alpha=0.6, color=colors[label], zorder=3)
            y_r = resid_model.params["Intercept"] + resid_model.params[pred] * x_p
            ax_resid.plot(x_p, y_r, color=colors[label], linewidth=2,
                          label=f"{label} (R²={r2_resid:.2f}, p={p_resid:.3f})")

        ax_raw.set_xlabel(pred_label, fontsize=14)
        ax_raw.set_ylabel(outcome_label, fontsize=14)
        ax_raw.set_title(f"{outcome_label} vs {pred_label}", fontsize=15)
        ax_raw.legend(fontsize=12)
        ax_raw.tick_params(labelsize=12)
        ax_raw.grid(True, alpha=0.4)

        ax_resid.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax_resid.set_xlabel(pred_label, fontsize=14)
        ax_resid.set_ylabel(
            f"{outcome_label} residuals\n(dist + impedance removed)", fontsize=14)
        ax_resid.set_title(
            f"{outcome_label} residuals ~ {pred_label}", fontsize=15)
        ax_resid.legend(fontsize=12)
        ax_resid.tick_params(labelsize=12)
        ax_resid.grid(True, alpha=0.4)

        fig.suptitle(f"{pred_label} — {outcome_label}", fontsize=16)
        plt.tight_layout()
        plt.show()
