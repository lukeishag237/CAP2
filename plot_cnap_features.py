"""
plot_cnap_features.py

Plots CNAP features from stored_df.pkl for SN and tibial myelinated fibers.
Uses Plotly for interactive plots with:
- Different shapes for each channel (circle, square, diamond)
- Different symbols for polarity (filled = positive, open = negative)
- Expdate labels on hover

For each outcome (peak_amplitude, peak_time), generates:
  Stage 1  — outcome ~ stim_rec_dist
  Stage 1b — outcome ~ impedance_kOhm
  Stage 2  — dist+imp residuals vs each histology predictor
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import statsmodels.formula.api as smf
import yaml
import os
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
CONFIG_FILE = "local_config.yml"
with open(CONFIG_FILE) as f:
    cfg = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = cfg["local_output_path"]
RAW_DATA_PATH     = cfg["raw_data_path"]

# ── Load stored_df.pkl ────────────────────────────────────────────────────────
df_all = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))
df_all["expdate"] = df_all["expdate"].astype(str)

# Peak time correction
df_all.loc[df_all["dataset_stim_polarity"] == 1, "peak_time"] += df_all["stim_PW"] / 2
df_all.loc[df_all["dataset_stim_polarity"] == 2, "peak_time"] -= df_all["stim_PW"] / 2


# Filter for SN myelinated fibers
sn_df = df_all[(df_all["rec_nerve"] == "SN") &
               (df_all["target_fiber_types"] == "myel") &
               (df_all["stim_rec_dist"] > 0) &
               (df_all["peak_amplitude"] > 0)].copy()

# Filter for tibial myelinated fibers
tibial_df = df_all[(df_all["rec_nerve"] == "tibial") &
                   (df_all["target_fiber_types"] == "myel") &
                   (df_all["stim_rec_dist"] > 0) &
                   (df_all["peak_amplitude"] > 0)].copy()

cvn_df = df_all[(df_all["rec_nerve"] == "cVN") &
               (df_all["target_fiber_types"] == "myel") &
               (df_all["stim_rec_dist"] > 0) &
               (df_all["peak_amplitude"] > 0)].copy()

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


sn_df = merge_all(apply_offsets(sn_df))
tibial_df = merge_all(apply_offsets(tibial_df))
cvn_df = merge_all(apply_offsets(cvn_df))

print(f"SN rows: {len(sn_df)}  |  Tibial rows: {len(tibial_df)}  |  cVN rows: {len(cvn_df)}")

# ── Plotting config ───────────────────────────────────────────────────────────
COLOR_SN     = "#2ca02c"  # Green
COLOR_TIBIAL = "#1f77b4"  # Blue
COLOR_CVN = "#FF0000" #Red

# Symbol mapping for channels (shapes)
CHANNEL_SYMBOLS = {
    1: "circle",
    2: "square", 
    3: "diamond"
}

# Symbol mapping for polarity (filled vs open)
# Positive polarity (1.0) = filled symbol
# Negative polarity (-1.0) = open symbol
def get_symbol_suffix(polarity):
    """Returns '-open' for negative polarity, '' for positive."""
    if polarity == -1.0:
        return "-open"
    return ""

DATASETS = [
    ("SN",     sn_df,     COLOR_SN,     "dash"),
    ("Tibial", tibial_df, COLOR_TIBIAL, "solid"),
    ("cVN", cvn_df, COLOR_CVN, 'dot'),
]

OUTCOMES = {"peak_amplitude": "Peak Amplitude (µV)",
            "peak_time":      "Peak Time (ms)"}

EXTRA_PREDICTORS = {
        "total_area_um2":           "Inner Area (µm²)",
        "hist_axon_count":          "Axon Count",
}

#     "total_area_um2":           "Inner Area (µm²)",
#     "area_ratio":               "Area Ratio (f/i)",
#     "avg_fascicle_size":        "Average Fascicle Size (µm²)",
#     "outer_thickness_um":       "Outer Thickness (µm)",
#     "area_proportion":          "Area Proportion (inner/outer)",
#     "distance_to_stim_cuff_um": "Electrode Fit (µm)",
# }


def add_scatter_traces(fig, df, outcome, pred, color, name, show_legend=True):
    """Add scatter traces for each channel and polarity combination."""
    df_v = df.dropna(subset=[outcome, pred])
    if len(df_v) < 1:
        return
    
    # Group by channel and polarity
    for ch in sorted(df_v['ch_num_neural'].unique()):
        for pol in sorted(df_v['dataset_stim_polarity'].unique()):
            df_subset = df_v[(df_v['ch_num_neural'] == ch) & 
                             (df_v['dataset_stim_polarity'] == pol)]
            
            if len(df_subset) == 0:
                continue
            
            # Get symbol
            base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
            symbol_suffix = get_symbol_suffix(pol)
            symbol = base_symbol + symbol_suffix
            
            # Create hover text with expdate
            hover_text = [f"Expdate: {row['expdate']}<br>"
                         f"Channel: {int(row['ch_num_neural'])}<br>"
                         f"Polarity: {'+' if row['dataset_stim_polarity'] == 1.0 else '-'}<br>"
                         f"{pred}: {row[pred]:.2f}<br>"
                         f"{outcome}: {row[outcome]:.2f}"
                         for _, row in df_subset.iterrows()]
            
            # Legend label (only show once per nerve type)
            legend_label = f"{name} (Ch{ch}, {'+'if pol==1.0 else '-'})"
            
            fig.add_trace(go.Scatter(
                x=df_subset[pred],
                y=df_subset[outcome],
                mode='markers',
                marker=dict(
                    symbol=symbol,
                    size=10,
                    color=color,
                    line=dict(width=1.5, color=color) if symbol_suffix == "-open" else dict(width=0)
                ),
                name=legend_label,
                legendgroup=name,
                showlegend=show_legend,
                hovertext=hover_text,
                hoverinfo='text'
            ))


def add_fit_line(fig, df, outcome, pred, color, name, linestyle="solid"):
    """Fit linear model and add regression line."""
    df_v = df.dropna(subset=[outcome, pred])
    if len(df_v) < 3:
        return
    
    model = smf.ols(f"{outcome} ~ {pred}", data=df_v).fit()
    r2 = model.rsquared
    p = model.pvalues.get(pred, np.nan)
    
    # Generate fit line
    xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
    yv = model.params["Intercept"] + model.params[pred] * xv
    
    fig.add_trace(go.Scatter(
        x=xv,
        y=yv,
        mode='lines',
        line=dict(color=color, width=2, dash=linestyle),
        name=f"{name} fit (R²={r2:.2f}, p={p:.3f})",
        legendgroup=name,
        showlegend=True,
        hoverinfo='skip'
    ))


def add_resid_scatter_traces(fig, df, outcome, pred, color, name, show_legend=True):
    """Add scatter traces for residuals after removing dist+imp effects."""
    df_v = df.dropna(subset=[outcome, "stim_rec_dist", "impedance_kOhm", pred])
    if len(df_v) < 3:
        return
    
    # Base model
    base = smf.ols(f"{outcome} ~ stim_rec_dist + impedance_kOhm", data=df_v).fit()
    df_v = df_v.copy()
    df_v["resid"] = base.resid
    
    # Group by channel and polarity
    for ch in sorted(df_v['ch_num_neural'].unique()):
        for pol in sorted(df_v['dataset_stim_polarity'].unique()):
            df_subset = df_v[(df_v['ch_num_neural'] == ch) & 
                             (df_v['dataset_stim_polarity'] == pol)]
            
            if len(df_subset) == 0:
                continue
            
            # Get symbol
            base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
            symbol_suffix = get_symbol_suffix(pol)
            symbol = base_symbol + symbol_suffix
            
            # Create hover text
            hover_text = [f"Expdate: {row['expdate']}<br>"
                         f"Channel: {int(row['ch_num_neural'])}<br>"
                         f"Polarity: {'+' if row['dataset_stim_polarity'] == 1.0 else '-'}<br>"
                         f"{pred}: {row[pred]:.2f}<br>"
                         f"Residual: {row['resid']:.2f}"
                         for _, row in df_subset.iterrows()]
            
            legend_label = f"{name} (Ch{ch}, {'+'if pol==1.0 else '-'})"
            
            fig.add_trace(go.Scatter(
                x=df_subset[pred],
                y=df_subset["resid"],
                mode='markers',
                marker=dict(
                    symbol=symbol,
                    size=10,
                    color=color,
                    line=dict(width=1.5, color=color) if symbol_suffix == "-open" else dict(width=0)
                ),
                name=legend_label,
                legendgroup=name,
                showlegend=show_legend,
                hovertext=hover_text,
                hoverinfo='text'
            ))


def add_resid_fit_line(fig, df, outcome, pred, color, name, linestyle="solid"):
    """Fit residual model and add regression line."""
    df_v = df.dropna(subset=[outcome, "stim_rec_dist", "impedance_kOhm", pred])
    if len(df_v) < 3:
        return
    
    # Base model
    base = smf.ols(f"{outcome} ~ stim_rec_dist + impedance_kOhm", data=df_v).fit()
    df_v = df_v.copy()
    df_v["resid"] = base.resid
    
    # Residual model
    rm = smf.ols(f"resid ~ {pred}", data=df_v).fit()
    r2 = rm.rsquared
    p = rm.pvalues.get(pred, np.nan)
    
    # Generate fit line
    xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
    yv = rm.params["Intercept"] + rm.params[pred] * xv
    
    fig.add_trace(go.Scatter(
        x=xv,
        y=yv,
        mode='lines',
        line=dict(color=color, width=2, dash=linestyle),
        name=f"{name} fit (R²={r2:.2f}, p={p:.3f})",
        legendgroup=name,
        showlegend=True,
        hoverinfo='skip'
    ))


# ── Stage 1: distance baseline ────────────────────────────────────────────────
print("\n=== Stage 1: Distance baseline ===")
for outcome, outcome_label in OUTCOMES.items():
    fig = go.Figure()
    
    for name, df, col, ls in DATASETS:
        add_scatter_traces(fig, df, outcome, "stim_rec_dist", col, name, show_legend=False)
        add_fit_line(fig, df, outcome, "stim_rec_dist", col, name, linestyle=ls)
    
    fig.update_layout(
        title=f"{outcome_label} ~ Distance",
        xaxis_title="Conduction Distance (mm)",
        yaxis_title=outcome_label,
        hovermode='closest',
        width=1000,
        height=700,
        template='plotly_white',
        font=dict(size=12)
    )
    
    fig.show()

# ── Stage 1b: impedance baseline ─────────────────────────────────────────────
print("\n=== Stage 1b: Impedance baseline ===")
for outcome, outcome_label in OUTCOMES.items():
    fig = go.Figure()
    
    for name, df, col, ls in DATASETS:
        add_scatter_traces(fig, df, outcome, "impedance_kOhm", col, name, show_legend=False)
        add_fit_line(fig, df, outcome, "impedance_kOhm", col, name, linestyle=ls)
    
    fig.update_layout(
        title=f"{outcome_label} ~ Impedance",
        xaxis_title="Impedance (kΩ)",
        yaxis_title=outcome_label,
        hovermode='closest',
        width=1000,
        height=700,
        template='plotly_white',
        font=dict(size=12)
    )
    
    fig.show()

# ── Stage 2: dist+imp residuals vs predictors ────────────────────────────────
print("\n=== Stage 2: Dist+Imp residuals ===")
for pred, pred_label in EXTRA_PREDICTORS.items():
    for outcome, outcome_label in OUTCOMES.items():
        # Create subplot with 2 columns
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=[f"{outcome_label} vs {pred_label}", 
                          f"Residuals ~ {pred_label}"],
            horizontal_spacing=0.12
        )
        
        # Left plot: raw data
        for name, df, col, ls in DATASETS:
            # Add scatter (only to first subplot)
            df_v = df.dropna(subset=[outcome, pred])
            if len(df_v) >= 1:
                for ch in sorted(df_v['ch_num_neural'].unique()):
                    for pol in sorted(df_v['dataset_stim_polarity'].unique()):
                        df_subset = df_v[(df_v['ch_num_neural'] == ch) & 
                                         (df_v['dataset_stim_polarity'] == pol)]
                        
                        if len(df_subset) == 0:
                            continue
                        
                        base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
                        symbol_suffix = get_symbol_suffix(pol)
                        symbol = base_symbol + symbol_suffix
                        
                        hover_text = [f"Expdate: {row['expdate']}<br>"
                                     f"Channel: {int(row['ch_num_neural'])}<br>"
                                     f"Polarity: {'+' if row['dataset_stim_polarity'] == 1.0 else '-'}<br>"
                                     f"{pred}: {row[pred]:.2f}<br>"
                                     f"{outcome}: {row[outcome]:.2f}"
                                     for _, row in df_subset.iterrows()]
                        
                        fig.add_trace(go.Scatter(
                            x=df_subset[pred],
                            y=df_subset[outcome],
                            mode='markers',
                            marker=dict(
                                symbol=symbol,
                                size=8,
                                color=col,
                                line=dict(width=1.5, color=col) if symbol_suffix == "-open" else dict(width=0)
                            ),
                            name=f"{name} (Ch{ch}, {'+'if pol==1.0 else '-'})",
                            legendgroup=name,
                            showlegend=False,
                            hovertext=hover_text,
                            hoverinfo='text'
                        ), row=1, col=1)
            
            # Add fit line
            if len(df_v) >= 3:
                model = smf.ols(f"{outcome} ~ {pred}", data=df_v).fit()
                r2 = model.rsquared
                p = model.pvalues.get(pred, np.nan)
                xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
                yv = model.params["Intercept"] + model.params[pred] * xv
                
                fig.add_trace(go.Scatter(
                    x=xv, y=yv,
                    mode='lines',
                    line=dict(color=col, width=2, dash=ls),
                    name=f"{name} fit (R²={r2:.2f}, p={p:.3f})",
                    legendgroup=name,
                    showlegend=True,
                    hoverinfo='skip'
                ), row=1, col=1)
        
        # Right plot: residuals
        for name, df, col, ls in DATASETS:
            df_v = df.dropna(subset=[outcome, "stim_rec_dist", "impedance_kOhm", pred])
            if len(df_v) < 3:
                continue
            
            base = smf.ols(f"{outcome} ~ stim_rec_dist + impedance_kOhm", data=df_v).fit()
            df_v = df_v.copy()
            df_v["resid"] = base.resid
            
            # Add scatter
            for ch in sorted(df_v['ch_num_neural'].unique()):
                for pol in sorted(df_v['dataset_stim_polarity'].unique()):
                    df_subset = df_v[(df_v['ch_num_neural'] == ch) & 
                                     (df_v['dataset_stim_polarity'] == pol)]
                    
                    if len(df_subset) == 0:
                        continue
                    
                    base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
                    symbol_suffix = get_symbol_suffix(pol)
                    symbol = base_symbol + symbol_suffix
                    
                    hover_text = [f"Expdate: {row['expdate']}<br>"
                                 f"Channel: {int(row['ch_num_neural'])}<br>"
                                 f"Polarity: {'+' if row['dataset_stim_polarity'] == 1.0 else '-'}<br>"
                                 f"{pred}: {row[pred]:.2f}<br>"
                                 f"Residual: {row['resid']:.2f}"
                                 for _, row in df_subset.iterrows()]
                    
                    fig.add_trace(go.Scatter(
                        x=df_subset[pred],
                        y=df_subset["resid"],
                        mode='markers',
                        marker=dict(
                            symbol=symbol,
                            size=8,
                            color=col,
                            line=dict(width=1.5, color=col) if symbol_suffix == "-open" else dict(width=0)
                        ),
                        name=f"{name} (Ch{ch}, {'+'if pol==1.0 else '-'})",
                        legendgroup=name,
                        showlegend=False,
                        hovertext=hover_text,
                        hoverinfo='text'
                    ), row=1, col=2)
            
            # Add fit line
            rm = smf.ols(f"resid ~ {pred}", data=df_v).fit()
            r2 = rm.rsquared
            p = rm.pvalues.get(pred, np.nan)
            xv = np.linspace(df_v[pred].min(), df_v[pred].max(), 200)
            yv = rm.params["Intercept"] + rm.params[pred] * xv
            
            fig.add_trace(go.Scatter(
                x=xv, y=yv,
                mode='lines',
                line=dict(color=col, width=2, dash=ls),
                name=f"{name} resid fit (R²={r2:.2f}, p={p:.3f})",
                legendgroup=name,
                showlegend=True,
                hoverinfo='skip'
            ), row=1, col=2)
        
        # Add horizontal line at y=0 for residuals
        fig.add_hline(y=0, line_dash="dash", line_color="gray", line_width=1, row=1, col=2)
        
        # Update axes
        fig.update_xaxes(title_text=pred_label, row=1, col=1)
        fig.update_xaxes(title_text=pred_label, row=1, col=2)
        fig.update_yaxes(title_text=outcome_label, row=1, col=1)
        fig.update_yaxes(title_text=f"{outcome_label} residuals<br>(dist+imp removed)", row=1, col=2)
        
        fig.update_layout(
            title_text=f"SN + Tibial — {pred_label} — {outcome_label}",
            hovermode='closest',
            width=1600,
            height=700,
            template='plotly_white',
            font=dict(size=11)
        )
        
        fig.show()

print("\n=== All plots complete ===")
