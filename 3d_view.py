"""
3d_view.py

3D scatter plots of CNAP features vs. conduction distance and impedance.

Figure 1: Peak Amplitude vs. Distance vs. Impedance
Figure 2: Peak Time      vs. Distance vs. Impedance

Data loading mirrors fit_distance_models.py exactly.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401 — registers 3D projection
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

cnap_df = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))

# ============================================================
# Load histology data
# ============================================================
hist_df = pd.read_csv(os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv"))

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

hist_summary = (
    hist_df.groupby(["expdate", "NerveCanonical"], as_index=False)
           .agg(hist_axon_count=("MinFeret", "count"))
)

# ============================================================
# Load fascicle data (needed for merge, even if not plotted)
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
    [["expdate", "NerveCanonical", "num_fascicles", "total_area_um2",
      "axon_area_um2", "area_ratio", "outer_area_um2", "outer_perimeter_um",
      "outer_thickness_um", "whole_area_um2"]]
    .copy()
)
fas_summary["area_proportion"] = (
    fas_summary["total_area_um2"] / fas_summary["outer_area_um2"]
)

# ============================================================
# Merge CNAP + histology + fascicle
# ============================================================
NERVES = ["cVN", "SN", "tibial"]

df = cnap_df[(cnap_df["rec_nerve"].isin(NERVES)) & (cnap_df["stim_rec_dist"] > 0)].copy()

_meta = pd.read_excel(os.path.join(RAW_DATA_PATH, "AnalysisScript_CAPs_Batch_parameters.xlsx"))
_meta["expdate"] = _meta["file_prefix_outputs_values"].str[:8]
_elec = (
    _meta.groupby(["expdate", "rec_nerve"], as_index=False)
         .agg(dataset_electrodetype=("dataset_electrodetype", "first"))
)
df = df.merge(_elec, on=["expdate", "rec_nerve"], how="left")

def channel_offset(row):
    n  = row["dataset_electrodetype"]
    ch = row["ch_num_neural"]
    if n == 2:
        return {1: -0.75, 2: +0.75}.get(ch, 0.0)
    elif n == 3:
        return {1: -1.5, 2: 0.0, 3: +1.5}.get(ch, 0.0)
    return 0.0

df["stim_rec_dist"] = df["stim_rec_dist"] + df.apply(channel_offset, axis=1)

df = pd.merge(df, hist_summary,
              left_on=["expdate", "rec_nerve"],
              right_on=["expdate", "NerveCanonical"], how="inner")
df = pd.merge(df, fas_summary, on=["expdate", "NerveCanonical"], how="inner")

# ============================================================
# Load impedance data
# ============================================================
imp_raw = pd.read_csv("all_impedences.csv")
imp_raw["expdate"] = imp_raw["expdate"].astype(str)

_hist_nerve = (
    hist_df[["Name", "NerveCanonical", "expdate"]]
    .drop_duplicates()
    .rename(columns={"Name": "Histology ID"})
)
imp_raw = imp_raw.merge(_hist_nerve, on=["Histology ID", "expdate"], how="left")

imp_long = pd.melt(
    imp_raw,
    id_vars=["expdate", "NerveCanonical"],
    value_vars=["ch_1_vs_ref", "ch_2_vs_ref", "ch_3_vs_ref"],
    var_name="ch_col", value_name="impedance_kOhm"
)
imp_long["ch_num_neural"] = imp_long["ch_col"].map(
    {"ch_1_vs_ref": 1, "ch_2_vs_ref": 2, "ch_3_vs_ref": 3}
)
imp_long = imp_long.drop(columns="ch_col").dropna(subset=["NerveCanonical"])

df = df.merge(
    imp_long[["expdate", "NerveCanonical", "ch_num_neural", "impedance_kOhm"]],
    on=["expdate", "NerveCanonical", "ch_num_neural"], how="left"
)

# ============================================================
# Prepare variables
# ============================================================
df = df[df["peak_amplitude"] > 0].copy()
df["fiber"]  = df["target_fiber_types"].astype("category")
df["nerve"]  = df["rec_nerve"].astype("category")

def classify_cvn(row):
    if row["nerve"] != "cVN":
        return row["nerve"]
    return "cVN_myel" if row["fiber"] == "myel" else "cVN_unmyel"

df["nerve_fiber_group"] = df.apply(classify_cvn, axis=1)

nerve_groups = {"tibial": "tibial", "SN": "SN", "cVN_myel": "cVN_myel"}
colors       = {"tibial": "tab:blue", "SN": "tab:green", "cVN_myel": "tab:red"}
markers      = {"tibial": "o",        "SN": "s",         "cVN_myel": "^"}

# ============================================================
# Exclusions — add (expdate, rec_nerve) tuples to remove
# specific experiments before fitting. Leave empty to use all.
# Examples:
#   ("20220531", "SN")
#   ("20221108", "tibial")
# ============================================================
EXCLUDE = [
    ("20220720", "tibial"),
    ("20221201", "SN")
]

if EXCLUDE:
    mask = pd.Series(
        [True] * len(df), index=df.index
    )
    for expdate, nerve in EXCLUDE:
        mask &= ~((df["expdate"] == expdate) & (df["rec_nerve"] == nerve))
    df = df[mask].copy()
    print(f"Excluded {(~mask).sum()} rows: {EXCLUDE}")

# ============================================================
# Fit distance + impedance plane for each nerve group
# ============================================================
def fit_plane(df_valid, outcome):
    """Fit outcome ~ dist + impedance and return a meshgrid surface."""
    model = smf.ols(f"{outcome} ~ stim_rec_dist + impedance_kOhm",
                    data=df_valid).fit()
    d_range = np.linspace(df_valid["stim_rec_dist"].min(),
                           df_valid["stim_rec_dist"].max(), 30)
    i_range = np.linspace(df_valid["impedance_kOhm"].min(),
                           df_valid["impedance_kOhm"].max(), 30)
    D, I = np.meshgrid(d_range, i_range)
    Z = (model.params["Intercept"]
         + model.params["stim_rec_dist"]  * D
         + model.params["impedance_kOhm"] * I)
    r2    = model.rsquared
    p_d   = model.pvalues.get("stim_rec_dist",  np.nan)
    p_i   = model.pvalues.get("impedance_kOhm", np.nan)
    return D, I, Z, r2, p_d, p_i

# ============================================================
# Plot
# ============================================================
OUTCOMES = {
    "peak_amplitude": "Peak Amplitude (µV)",
    "peak_time":      "Peak Time (ms)",
    "auc":            "AUC (µV·ms)",
}

for outcome, outcome_label in OUTCOMES.items():
    fig = plt.figure(figsize=(11, 8))
    ax  = fig.add_subplot(111, projection="3d")

    for label, code in nerve_groups.items():
        df_sub   = df[df["nerve_fiber_group"] == code].copy()
        df_valid = df_sub.dropna(subset=[outcome, "stim_rec_dist", "impedance_kOhm"])
        if len(df_valid) < 3:
            continue

        # Scatter points
        ax.scatter(
            df_valid["stim_rec_dist"],
            df_valid["impedance_kOhm"],
            df_valid[outcome],
            color=colors[label],
            marker=markers[label],
            s=50, alpha=0.8,
            label=label, zorder=5
        )

        # Fitted plane
        D, I, Z, r2, p_d, p_i = fit_plane(df_valid, outcome)
        ax.plot_surface(D, I, Z,
                        color=colors[label], alpha=0.18, linewidth=0)
        # Invisible proxy for the legend entry with stats
        ax.plot([], [], [], color=colors[label], linewidth=3,
                label=f"{label}  R²={r2:.2f}  p(dist)={p_d:.3f}  p(imp)={p_i:.3f}")

    ax.set_xlabel("Conduction Distance (mm)", fontsize=12, labelpad=10)
    ax.set_ylabel("Impedance (kΩ)",           fontsize=12, labelpad=10)
    ax.set_zlabel(outcome_label,              fontsize=12, labelpad=10)
    ax.set_title(f"{outcome_label}\nvs Distance & Impedance", fontsize=14)

    # Deduplicate legend (scatter + plane proxy both add entries)
    handles, labels_ = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels_):
        if l not in seen:
            seen[l] = h
    ax.legend(seen.values(), seen.keys(), fontsize=9, loc="upper left")

    plt.tight_layout()
    plt.show()
