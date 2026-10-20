"""
compare_artifact_removal.py

Compares two datasets for SN myelinated fibers:
  - "Original"  : stored_df.pkl  (pipeline without artifact removal)
  - "Artefact-removed": stored_df_artifact_removed_SN.mat

For each, fits:
  Stage 1  — peak_amplitude / peak_time ~ stim_rec_dist
  Stage 1b — peak_amplitude / peak_time ~ impedance_kOhm
  Stage 2  — dist+imp residuals vs each histology predictor

All figures show Original (solid) vs Artifact-removed (dashed) side by side.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import statsmodels.formula.api as smf
import scipy.io as sio
import yaml
import os
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
CONFIG_FILE = "local_config.yml"
with open(CONFIG_FILE) as f:
    cfg = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = cfg["local_output_path"]
RAW_DATA_PATH     = cfg["raw_data_path"]

# ── Load original stored_df.pkl ───────────────────────────────────────────────
orig_df = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))
orig_df["expdate"] = orig_df["expdate"].astype(str)
orig_df = orig_df[(orig_df["rec_nerve"] == "SN") &
                  (orig_df["target_fiber_types"] == "myel") &
                  (orig_df["stim_rec_dist"] > 0) &
                  (orig_df["peak_amplitude"] > 0)].copy()

# ── Load original tibial from stored_df.pkl ───────────────────────────────────
orig_tib_df = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))
orig_tib_df["expdate"] = orig_tib_df["expdate"].astype(str)

orig_tib_df = orig_tib_df[(orig_tib_df["rec_nerve"] == "tibial") &
                          (orig_tib_df["target_fiber_types"] == "myel") &
                          (orig_tib_df["stim_rec_dist"] > 0) &
                          (orig_tib_df["peak_amplitude"] > 0)].copy()


# ── Load artifact-removed mat ─────────────────────────────────────────────────
mat = sio.loadmat(
    os.path.join(LOCAL_OUTPUT_PATH, "stored_df_artifact_removed_SN.mat"),
    squeeze_me=True, struct_as_record=False)
rows = mat["all_rows"]

art_records = []
for r in rows:
    art_records.append({
        "expdate":               str(int(r.expdate)),
        "rec_nerve":             str(r.rec_nerve),
        "target_fiber_types":    str(r.target_fiber_types),
        "ch_num_neural":         int(r.ch_num_neural),
        "stim_rec_dist":         float(r.stim_rec_dist),
        "peak_amplitude":        float(r.peak_amplitude),
        "peak_time":             float(r.peak_time),
        "auc":                   float(r.auc),
        "dataset_stim_polarity": float(r.dataset_stim_polarity),
    })
art_df = pd.DataFrame(art_records)
art_df = art_df[art_df["peak_amplitude"] > 0].copy()

# ── Load tibial artifact-removed mat ─────────────────────────────────────────
mat_tib = sio.loadmat(
    os.path.join(LOCAL_OUTPUT_PATH, "stored_df_artifact_removed_tibial.mat"),
    squeeze_me=True, struct_as_record=False)
rows_tib = mat_tib["all_rows"]

tib_records = []
for r in rows_tib:
    tib_records.append({
        "expdate":               str(int(r.expdate)),
        "rec_nerve":             str(r.rec_nerve),
        "target_fiber_types":    str(r.target_fiber_types),
        "ch_num_neural":         int(r.ch_num_neural),
        "stim_rec_dist":         float(r.stim_rec_dist),
        "peak_amplitude":        float(r.peak_amplitude),
        "peak_time":             float(r.peak_time),
        "auc":                   float(r.auc),
        "dataset_stim_polarity": float(r.dataset_stim_polarity),
    })

tib_df = pd.DataFrame(tib_records)
tib_df = tib_df[tib_df["peak_amplitude"] > 0].copy()


# ── Shared histology + fascicle + impedance loading ───────────────────────────
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

hist_df = pd.read_csv(os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv"))
hist_df["expdate"]        = hist_df["Name"].map(mapping)
hist_df["NerveCanonical"] = hist_df["NerveType"].map(nerve_map)
hist_df = hist_df.dropna(subset=["expdate","NerveCanonical"])
hist_df["expdate"] = hist_df["expdate"].astype(str)

hist_summary = (hist_df.groupby(["expdate","NerveCanonical"], as_index=False)
                .agg(hist_axon_count=("MinFeret","count"),
                     hist_mean_diam=("MinFeret","mean")))

BASE_SUMMARY_DIR = os.path.join(RAW_DATA_PATH, "..", "60x", "derivatives", "summary")
fas_df = pd.read_csv(os.path.join(BASE_SUMMARY_DIR, "fascicle_data.csv"))
fas_df["expdate"] = fas_df["subject"].map(mapping)
fas_df["expdate"] = fas_df["expdate"].astype(str)
nlu = hist_df[["Name","NerveType"]].drop_duplicates()
fas_df = fas_df.merge(nlu, left_on="subject", right_on="Name", how="left")
fas_df["NerveCanonical"] = fas_df["NerveType"].map(nerve_map)
fas_summary = (
    fas_df[fas_df["mask_type"]=="summary"]
    [["expdate","NerveCanonical","num_fascicles","total_area_um2","area_ratio",
      "outer_area_um2","outer_perimeter_um","outer_thickness_um","whole_area_um2"]]
    .copy()
)
fas_summary["area_proportion"] = fas_summary["total_area_um2"] / fas_summary["outer_area_um2"]
fas_summary["avg_fascicle_size"] = fas_summary["total_area_um2"] / fas_summary["num_fascicles"]

dist_summary = (hist_df.groupby(["Name","expdate"], as_index=False)
                .agg(mean_dist=("distance_to_stim_cuff_um","mean")))

# Impedance
imp_raw = pd.read_csv("all_impedences.csv")
imp_raw["expdate"] = imp_raw["expdate"].astype(str)
_hn = (hist_df[["Name","NerveCanonical","expdate"]].drop_duplicates()
       .rename(columns={"Name":"Histology ID"}))
imp_raw = imp_raw.merge(_hn, on=["Histology ID","expdate"], how="left")
imp_long = pd.melt(imp_raw, id_vars=["expdate","NerveCanonical"],
                   value_vars=["ch_1_vs_ref","ch_2_vs_ref","ch_3_vs_ref"],
                   var_name="ch_col", value_name="impedance_kOhm")
imp_long["ch_num_neural"] = imp_long["ch_col"].map(
    {"ch_1_vs_ref":1,"ch_2_vs_ref":2,"ch_3_vs_ref":3})
imp_long = imp_long.drop(columns="ch_col").dropna(subset=["NerveCanonical"])

# ── Electrode type + channel offset ──────────────────────────────────────────
_meta = pd.read_excel(os.path.join(RAW_DATA_PATH,"AnalysisScript_CAPs_Batch_parameters.xlsx"))
_meta["expdate"] = _meta["file_prefix_outputs_values"].str[:8]
_elec = (_meta.groupby(["expdate","rec_nerve"], as_index=False)
         .agg(dataset_electrodetype=("dataset_electrodetype","first")))

def apply_offsets(df):
    df = df.merge(_elec, on=["expdate","rec_nerve"], how="left")
    def off(row):
        n, ch = row["dataset_electrodetype"], row["ch_num_neural"]
        return {2:{1:-0.75,2:+0.75}, 3:{1:-1.5,2:0.,3:+1.5}}.get(n,{}).get(ch,0.)
    df["stim_rec_dist"] = df["stim_rec_dist"] + df.apply(off, axis=1)
    return df

def merge_all(df):
    df = pd.merge(df, hist_summary,
                  left_on=["expdate","rec_nerve"], right_on=["expdate","NerveCanonical"], how="inner")
    df = pd.merge(df, fas_summary, on=["expdate","NerveCanonical"], how="inner")
    df = df.merge(dist_summary[["expdate","mean_dist"]], on="expdate", how="left")
    df = df.rename(columns={"mean_dist":"distance_to_stim_cuff_um"})
    df = df.merge(imp_long[["expdate","NerveCanonical","ch_num_neural","impedance_kOhm"]],
                  on=["expdate","NerveCanonical","ch_num_neural"], how="left")
    return df

# orig_df = merge_all(apply_offsets(orig_df))
# art_df  = merge_all(apply_offsets(art_df))
orig_df = merge_all(apply_offsets(orig_df)) # no apply_offsets
orig_tib_df = merge_all(apply_offsets(orig_tib_df))
art_df  = merge_all(art_df)
tib_df = merge_all(tib_df)



print(f"Original SN rows: {len(orig_df)}  |  Artifact-removed SN rows: {len(art_df)}")

# ── Plotting config ───────────────────────────────────────────────────────────
COLOR_ORIG = "tab:green"
COLOR_ART  = "tab:orange"
COLOR_TIB  = "tab:blue"
COLOR_ORIG_TIB = "tab:red"

DATASETS = [
    ("Original SN",             orig_df,      COLOR_ORIG,     "-"),
    ("Artifact-removed SN",     art_df,       COLOR_ART,      "--"),
    ("Original tibial",         orig_tib_df,  COLOR_ORIG_TIB, "-."),
    ("Artifact-removed tibial", tib_df,       COLOR_TIB,      ":")
]



OUTCOMES = {"peak_amplitude": "Peak Amplitude (µV)",
            "peak_time":      "Peak Time (ms)"}

EXTRA_PREDICTORS = {
    "hist_axon_count":          "Axon Count",
    "total_area_um2":           "Inner Area (µm²)",
    "area_ratio":               "Area Ratio (f/i)",
    "avg_fascicle_size":        "Average Fascicle Size (µm²)",
    "outer_thickness_um":       "Outer Thickness (µm)",
    "area_proportion":          "Area Proportion (inner/outer)",
    "distance_to_stim_cuff_um": "Electrode Fit (µm)",
}

legend_patches = [
    mlines.Line2D([], [], color=COLOR_ORIG,     linestyle="-",  linewidth=2, label="Original SN"),
    mlines.Line2D([], [], color=COLOR_ART,      linestyle="--", linewidth=2, label="Artifact-removed SN"),
    mlines.Line2D([], [], color=COLOR_ORIG_TIB, linestyle="-.", linewidth=2, label="Original tibial"),
    mlines.Line2D([], [], color=COLOR_TIB,      linestyle=":",  linewidth=2, label="Artifact-removed tibial"),
]



def fit_and_plot(ax, df, outcome, pred, color, ls):
    df_v = df.dropna(subset=[outcome, pred])
    if len(df_v) < 3:
        return
    model = smf.ols(f"{outcome} ~ {pred}", data=df_v).fit()
    r2 = model.rsquared
    p  = model.pvalues.get(pred, np.nan)
    ax.scatter(df_v[pred], df_v[outcome], color=color, alpha=0.45, s=25, zorder=3)
    xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
    ax.plot(xv, model.params["Intercept"] + model.params[pred]*xv,
            color=color, linestyle=ls, linewidth=2,
            label=f"R²={r2:.2f} p={p:.3f}")

def fit_resid_and_plot(ax, df, outcome, pred, color, ls):
    df_v = df.dropna(subset=[outcome, "stim_rec_dist", "impedance_kOhm", pred])
    if len(df_v) < 3:
        return
    base  = smf.ols(f"{outcome} ~ stim_rec_dist + impedance_kOhm", data=df_v).fit()
    df_v  = df_v.copy(); df_v["resid"] = base.resid
    rm    = smf.ols(f"resid ~ {pred}", data=df_v).fit()
    r2    = rm.rsquared; p = rm.pvalues.get(pred, np.nan)
    ax.scatter(df_v[pred], df_v["resid"], color=color, alpha=0.45, s=25, zorder=3)
    xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
    ax.plot(xv, rm.params["Intercept"] + rm.params[pred]*xv,
            color=color, linestyle=ls, linewidth=2,
            label=f"R²={r2:.2f} p={p:.3f}")

# ── Stage 1: distance baseline ────────────────────────────────────────────────
print("\n=== Stage 1: Distance baseline ===")
for outcome, outcome_label in OUTCOMES.items():
    fig, ax = plt.subplots(figsize=(8,6))
    for name, df, col, ls in DATASETS:
        fit_and_plot(ax, df, outcome, "stim_rec_dist", col, ls)
    ax.set_xlabel("Conduction Distance (mm)", fontsize=14)
    ax.set_ylabel(outcome_label, fontsize=14)
    ax.set_title(f"SN — {outcome_label} ~ Distance", fontsize=15)
    ax.legend(handles=legend_patches + ax.get_legend_handles_labels()[0], fontsize=11)
    ax.tick_params(labelsize=12); ax.grid(True, alpha=0.4)
    plt.tight_layout(); plt.show()

# ── Stage 1b: impedance baseline ─────────────────────────────────────────────
print("\n=== Stage 1b: Impedance baseline ===")
for outcome, outcome_label in OUTCOMES.items():
    fig, ax = plt.subplots(figsize=(8,6))
    for name, df, col, ls in DATASETS:
        fit_and_plot(ax, df, outcome, "impedance_kOhm", col, ls)
    ax.set_xlabel("Impedance (kΩ)", fontsize=14)
    ax.set_ylabel(outcome_label, fontsize=14)
    ax.set_title(f"SN — {outcome_label} ~ Impedance", fontsize=15)
    ax.legend(handles=legend_patches + ax.get_legend_handles_labels()[0], fontsize=11)
    ax.tick_params(labelsize=12); ax.grid(True, alpha=0.4)
    plt.tight_layout(); plt.show()

# ── Stage 2: dist+imp residuals vs predictors ────────────────────────────────
print("\n=== Stage 2: Dist+Imp residuals ===")
for pred, pred_label in EXTRA_PREDICTORS.items():
    for outcome, outcome_label in OUTCOMES.items():
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        ax_raw, ax_resid = axes

        for name, df, col, ls in DATASETS:
            fit_and_plot(ax_raw,   df, outcome, pred, col, ls)
            fit_resid_and_plot(ax_resid, df, outcome, pred, col, ls)

        ax_raw.set_xlabel(pred_label, fontsize=13)
        ax_raw.set_ylabel(outcome_label, fontsize=13)
        ax_raw.set_title(f"{outcome_label} vs {pred_label}", fontsize=13)
        ax_raw.legend(handles=legend_patches + ax_raw.get_legend_handles_labels()[0], fontsize=9)
        ax_raw.tick_params(labelsize=11); ax_raw.grid(True, alpha=0.4)

        ax_resid.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax_resid.set_xlabel(pred_label, fontsize=13)
        ax_resid.set_ylabel(f"{outcome_label} residuals\n(dist+imp removed)", fontsize=13)
        ax_resid.set_title(f"Residuals ~ {pred_label}", fontsize=13)
        ax_resid.legend(handles=legend_patches + ax_resid.get_legend_handles_labels()[0], fontsize=9)
        ax_resid.tick_params(labelsize=11); ax_resid.grid(True, alpha=0.4)

        fig.suptitle(f"SN + Tibial — {pred_label} — {outcome_label}", fontsize=15)
        plt.tight_layout(); plt.show()
