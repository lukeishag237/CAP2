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
HIST_CSV = os.path.join(
    RAW_DATA_PATH, "..", "60x", "derivatives", "summary",
    "all_samples_fiber_diameters.csv"
)
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
# Merge CNAP + histology
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

# ============================================================
# Prepare variables
# ============================================================
METRICS = {
    "auc": "log_auc",
    "peak_amplitude": "log_peak_amp",
    "peak_time": "log_peak_time"
}

for raw, log_name in METRICS.items():
    df = df[df[raw] > 0]
    df[log_name] = np.log(df[raw])

df["log_dist"] = np.log(df["stim_rec_dist"])
df["fiber"] = df["target_fiber_types"].astype("category")
df["nerve"] = df["rec_nerve"].astype("category")
df["subject"] = df["expdate"].astype("category")

# ============================================================
# NEW: Split cVN into myelinated vs unmyelinated
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
# Helper: mixed model log(y) ~ log(distance) + (1|subject)
# ============================================================
def fit_distance_mixed(df_sub, log_metric):
    model = smf.mixedlm(
        f"{log_metric} ~ log_dist",
        df_sub,
        groups=df_sub["subject"]
    )
    result = model.fit(method="lbfgs")

    fe_pred = result.predict(exog=df_sub[["log_dist"]])
    df_sub[f"{log_metric}_dist_me_resid"] = df_sub[log_metric] - fe_pred

    return result, df_sub

# ============================================================
# 1. Fit model separately for each nerve-fiber group
# ============================================================
nerve_groups = {
    "tibial": "tibial",
    "sciatic": "SN",
    "cVN_myel": "cVN_myel",
    "cVN_unmyel": "cVN_unmyel"
}

for metric, log_metric in METRICS.items():
    print(f"\n================ METRIC: {metric.upper()} =================\n")

    for label, code in nerve_groups.items():
        print(f"\n---- GROUP: {label.upper()} ----")

        df_sub = df[df["nerve_fiber_group"] == code].copy()
        if df_sub.empty:
            print("No data for this group.")
            continue

        result, df_sub = fit_distance_mixed(df_sub, log_metric)
        print(result.summary())

        plt.figure(figsize=(7,6))
        plt.scatter(
            df_sub["stim_rec_dist"],
            df_sub[f"{log_metric}_dist_me_resid"],
            alpha=0.6,
            label=label
        )
        plt.axhline(0, color="black", linestyle="--")
        plt.xlabel("Conduction distance (mm)")
        plt.ylabel(f"{metric} residual (distance + subject removed)")
        plt.title(f"{metric.upper()} — {label} — mixed model residuals")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

# ============================================================
# 2. Combined model: log(y) ~ log(distance) + nerve_fiber_group + (1|subject)
# ============================================================
for metric, log_metric in METRICS.items():
    print(f"\n================ COMBINED MODEL: {metric.upper()} =================\n")

    model = smf.mixedlm(
        f"{log_metric} ~ log_dist + nerve_fiber_group",
        df,
        groups=df["subject"]
    )
    result = model.fit(method="lbfgs")

    fe_pred = result.predict(exog=df[["log_dist", "nerve_fiber_group"]])
    df[f"{log_metric}_combined_me_resid"] = df[log_metric] - fe_pred

    print(result.summary())

    plt.figure(figsize=(7,6))
    for grp in df["nerve_fiber_group"].unique():
        dfn = df[df["nerve_fiber_group"] == grp]
        plt.scatter(
            dfn["stim_rec_dist"],
            dfn[f"{log_metric}_combined_me_resid"],
            alpha=0.6,
            label=grp
        )
    plt.axhline(0, color="black", linestyle="--")
    plt.xlabel("Conduction distance (mm)")
    plt.ylabel(f"{metric} residual (distance + group + subject removed)")
    plt.title(f"{metric.upper()} — combined mixed model residuals")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# ============================================================
# 3. Second-stage model: how histology explains remaining variance
# ============================================================
for metric, log_metric in METRICS.items():
    print(f"\n================ SECOND-STAGE HISTOLOGY MODEL: {metric.upper()} =================\n")

    model2 = smf.ols(
        f"{log_metric}_combined_me_resid ~ hist_mean_diam + hist_median_diam + hist_total_diam + hist_axon_count",
        data=df
    ).fit()

    print(model2.summary())

    plt.figure(figsize=(7,6))
    plt.scatter(df["hist_axon_count"], df[f"{log_metric}_combined_me_resid"], alpha=0.6)
    plt.axhline(0, color="black", linestyle="--")
    plt.xlabel("Axon Count")
    plt.ylabel(f"{metric} residual")
    plt.title(f"{metric.upper()} residuals vs histology (mean diameter)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
