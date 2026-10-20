"""
predict_distance.py

Predict conduction distance from peak amplitude and peak time.
Fits linear models: conduction_distance ~ peak_amplitude (or peak_time)
for SN, tibial, and cVN separately.

Accounts for electrode offset based on electrode type and channel:
- Bipolar (ch 1,2): ch1 -= 0.75mm, ch2 += 0.75mm
- Tripolar (ch 1,2,3): ch1 -= 1.5mm, ch2 += 0mm, ch3 += 1.5mm

Plots predictions vs actual distances, then residuals vs predictor.
Uses Plotly with different symbols for polarity (+ for positive, dot for negative).
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import statsmodels.formula.api as smf
import yaml
import os
from pathlib import Path
import sys

# Add parent directory to path to import cnap_analyzer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnap_analyzer import CNAPDataProcessor

# ── Config ────────────────────────────────────────────────────────────────────
CONFIG_FILE = "local_config.yml"
with open(CONFIG_FILE) as f:
    cfg = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = cfg["local_output_path"]
RAW_DATA_PATH     = cfg["raw_data_path"]

# ── Load stored_df.pkl ────────────────────────────────────────────────────────
df_all = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))
df_all["expdate"] = df_all["expdate"].astype(str)

# If max_stim_amplitude_mA is missing, delete the pickle and regenerate
if "max_stim_amplitude_mA" not in df_all.columns:
    print("Warning: max_stim_amplitude_mA not found in stored_df.pkl")
    print("Regenerating stored_df.pkl from preprocessed data...")
    
    import glob
    from pathlib import Path
    
    # Reload data from cnap_analyzer
    processor = CNAPDataProcessor(str(Path(RAW_DATA_PATH) / 'AnalysisScript_CAPs_Batch_parameters.xlsx'), LOCAL_OUTPUT_PATH)
    
    # Use both paths for data loading
    ar_path = cfg.get("AR_preprocessed_CAP_path")
    processor.process_data(cfg.get("preprocessed_CAP_path"), ar_path, flag_plot=False)
    
    df_all = processor.data_df
    df_all.to_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))
    print("Regenerated stored_df.pkl with max_stim_amplitude_mA field")
else:
    print(f"Loaded stored_df.pkl with {len(df_all)} records")

# Import CNAPDataProcessor if we haven't already
try:
    processor
except NameError:
    pass  # Already imported above if needed

# Filter for myelinated fibers only
df_all = df_all[(df_all["target_fiber_types"] == "myel") &
                (df_all["stim_rec_dist"] > 0) &
                (df_all["peak_amplitude"] > 0) &
                (df_all["peak_time"] > 0)].copy()

# ── Load electrode type metadata ──────────────────────────────────────────────
_meta = pd.read_excel(os.path.join(RAW_DATA_PATH,"AnalysisScript_CAPs_Batch_parameters.xlsx"))
_meta["expdate"] = _meta["file_prefix_outputs_values"].str[:8]
_meta["rec_nerve"] = _meta["rec_nerve"]
_elec = (_meta.groupby(["expdate","rec_nerve"], as_index=False)
         .agg(dataset_electrodetype=("dataset_electrodetype","first")))

# Merge electrode type into main dataframe
df_all = df_all.merge(_elec, on=["expdate","rec_nerve"], how="left")

# ── Apply electrode offset corrections ────────────────────────────────────────
def apply_electrode_offset(row):
    """Apply offset based on electrode type and channel."""
    elec_type = row["dataset_electrodetype"]
    ch = row["ch_num_neural"]
    
    # Bipolar (type 2)
    if elec_type == 2:
        if ch == 1:
            return -0.75
        elif ch == 2:
            return 0.75
    # Tripolar (type 3)
    elif elec_type == 3:
        if ch == 1:
            return -1.5
        elif ch == 2:
            return 0.0
        elif ch == 3:
            return 1.5
    
    return 0.0

df_all["offset_mm"] = df_all.apply(apply_electrode_offset, axis=1)
df_all["stim_rec_dist_corrected"] = df_all["stim_rec_dist"] + df_all["offset_mm"]

# ── Symbol mapping for polarity ──────────────────────────────────────────────
def get_symbol_suffix(polarity):
    """Returns '-open' for negative polarity (-1), '' for positive (1)."""
    if polarity == -1.0:
        return "-open"
    return ""

CHANNEL_SYMBOLS = {
    1: "circle",
    2: "square", 
    3: "diamond"
}

# ── Plotting config ───────────────────────────────────────────────────────────
nerve_colors = {
    "SN": "#2ca02c",      # Green
    "tibial": "#1f77b4",  # Blue
    "cVN": "#d62728"      # Red
}

OUTCOMES = {
    "peak_amplitude": "Peak Amplitude (µV)",
    "peak_time": "Peak Time (ms)"
}

print(f"Total myelinated records: {len(df_all)}")
print(f"Records per nerve:")
for nerve in ["SN", "tibial", "cVN"]:
    count = len(df_all[df_all["rec_nerve"] == nerve])
    print(f"  {nerve}: {count}")


def fit_and_plot_prediction(predictor, predictor_label):
    """Fit models and plot predictions for each nerve."""
    print(f"\n=== Fitting conduction distance ~ {predictor_label} ===")
    
    fig = go.Figure()
    
    for nerve in ["SN", "tibial", "cVN"]:
        df_nerve = df_all[df_all["rec_nerve"] == nerve].dropna(subset=[predictor])
        
        if len(df_nerve) < 3:
            print(f"  {nerve}: Not enough data (n={len(df_nerve)})")
            continue
        
        # Fit linear model
        model = smf.ols(f"stim_rec_dist_corrected ~ {predictor}", data=df_nerve).fit()
        r2 = model.rsquared
        p = model.pvalues.get(predictor, np.nan)
        
        print(f"  {nerve}: n={len(df_nerve)}, R²={r2:.3f}, p={p:.4f}")
        
        # Store model in dataframe for residual calculation
        df_nerve["predicted_dist"] = model.predict(df_nerve)
        df_nerve["residual"] = df_nerve["stim_rec_dist_corrected"] - df_nerve["predicted_dist"]
        
        # Add actual points grouped by channel and polarity
        for ch in sorted(df_nerve['ch_num_neural'].unique()):
            for pol in sorted(df_nerve['dataset_stim_polarity'].unique()):
                df_subset = df_nerve[(df_nerve['ch_num_neural'] == ch) & 
                                     (df_nerve['dataset_stim_polarity'] == pol)]
                
                if len(df_subset) == 0:
                    continue
                
                base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
                symbol_suffix = get_symbol_suffix(pol)
                symbol = base_symbol + symbol_suffix
                
                hover_text = [f"Expdate: {row['expdate']}<br>"
                             f"Channel: {int(row['ch_num_neural'])}<br>"
                             f"Polarity: {'+' if row['dataset_stim_polarity'] == 1.0 else '-'}<br>"
                             f"{predictor}: {row[predictor]:.2f}<br>"
                             f"Distance: {row['stim_rec_dist_corrected']:.2f} mm"
                             for _, row in df_subset.iterrows()]
                
                fig.add_trace(go.Scatter(
                    x=df_subset[predictor],
                    y=df_subset["stim_rec_dist_corrected"],
                    mode='markers',
                    marker=dict(
                        symbol=symbol,
                        size=10,
                        color=nerve_colors[nerve],
                        line=dict(width=1.5, color=nerve_colors[nerve]) if symbol_suffix == "-open" else dict(width=0)
                    ),
                    name=f"{nerve} (Ch{ch}, {'+'if pol==1.0 else '-'})",
                    legendgroup=nerve,
                    showlegend=bool(ch == 1 and pol == 1.0),  # Convert to Python bool
                    hovertext=hover_text,
                    hoverinfo='text'
                ))
        
        # Add fit line
        xv = np.linspace(df_nerve[predictor].min(), df_nerve[predictor].max(), 200)
        yv = model.params["Intercept"] + model.params[predictor] * xv
        
        fig.add_trace(go.Scatter(
            x=xv,
            y=yv,
            mode='lines',
            line=dict(color=nerve_colors[nerve], width=2.5, dash="solid"),
            name=f"{nerve} fit (R²={r2:.2f})",
            legendgroup=nerve,
            showlegend=True,
            hoverinfo='skip'
        ))
    
    fig.update_layout(
        title=f"Conduction Distance ~ {predictor_label}",
        xaxis_title=predictor_label,
        yaxis_title="Conduction Distance (mm)",
        hovermode='closest',
        width=1100,
        height=700,
        template='plotly_white',
        font=dict(size=12)
    )
    
    fig.show()


def plot_residual_differences(predictor, predictor_label):
    """Plot residual differences between paired polarities (polarity 1 vs 2) vs stimulus amplitude.
    
    Stimulus amplitude is the MAXIMUM amplitude value from the list of stimulus amplitudes.
    Loaded from max_stim_amplitude_mA in the stored dataframe.
    """
    print(f"\n=== Plotting residual differences for {predictor_label} ===")
    
    fig = go.Figure()
    
    for nerve in ["SN", "tibial", "cVN"]:
        df_nerve = df_all[df_all["rec_nerve"] == nerve].dropna(subset=[predictor, "max_stim_amplitude_mA"])
        
        if len(df_nerve) < 3:
            continue
        
        # Fit model
        model = smf.ols(f"stim_rec_dist_corrected ~ {predictor}", data=df_nerve).fit()
        df_nerve = df_nerve.copy()
        df_nerve["predicted_dist"] = model.predict(df_nerve)
        df_nerve["residual"] = df_nerve["stim_rec_dist_corrected"] - df_nerve["predicted_dist"]
        
        # Create pair identifier (expdate, nerve, channel)
        df_nerve["pair_id"] = (df_nerve["expdate"].astype(str) + "_" + 
                               df_nerve["rec_nerve"].astype(str) + "_ch" + 
                               df_nerve["ch_num_neural"].astype(str))
        
        # Pivot to get polarity 1 and 2 side by side
        pivot = df_nerve.pivot_table(
            index="pair_id",
            columns="dataset_stim_polarity",
            values=["residual", "max_stim_amplitude_mA"],
            aggfunc="first"
        )
        
        # Flatten MultiIndex
        pivot.columns = [f"{metric}_pol{int(pol)}" if pol != 0 else metric 
                        for metric, pol in pivot.columns]
        
        # Drop incomplete pairs
        pivot = pivot.dropna()
        
        if len(pivot) == 0:
            print(f"  {nerve}: No complete polarity pairs")
            continue
        
        print(f"  {nerve}: {len(pivot)} complete pairs")
        
        # Calculate difference in residuals (pol 2 - pol 1)
        pivot["residual_diff"] = pivot["residual_pol2"] - pivot["residual_pol1"]
        
        # Use maximum stimulus amplitude (should be same for both polarities)
        pivot["stim_amplitude_mA"] = pivot["max_stim_amplitude_mA_pol1"]
        
        # Reset index to get pair info back
        pivot_reset = pivot.reset_index()
        
        # Extract channel info from pair_id for coloring
        pivot_reset["ch_num"] = pivot_reset["pair_id"].str.extract(r"ch(\d)")[0].astype(int)
        
        # Add points grouped by channel
        for ch in sorted(pivot_reset['ch_num'].unique()):
            df_subset = pivot_reset[pivot_reset['ch_num'] == ch]
            
            if len(df_subset) == 0:
                continue
            
            base_symbol = CHANNEL_SYMBOLS.get(ch, "circle")
            
            hover_text = [f"Pair: {row['pair_id']}<br>"
                         f"Channel: {int(row['ch_num'])}<br>"
                         f"Max Stimulus Amplitude: {row['stim_amplitude_mA']:.3f} mA<br>"
                         f"Residual Diff (Pol2-Pol1): {row['residual_diff']:.2f} mm<br>"
                         f"Pol1 residual: {row['residual_pol1']:.2f} mm<br>"
                         f"Pol2 residual: {row['residual_pol2']:.2f} mm"
                         for _, row in df_subset.iterrows()]
            
            fig.add_trace(go.Scatter(
                x=df_subset['stim_amplitude_mA'],
                y=df_subset['residual_diff'],
                mode='markers',
                marker=dict(
                    symbol=base_symbol,
                    size=10,
                    color=nerve_colors[nerve],
                    line=dict(width=1)
                ),
                name=f"{nerve} Ch{ch}",
                legendgroup=nerve,
                showlegend=True,
                hovertext=hover_text,
                hoverinfo='text'
            ))
    
    # Add zero line
    fig.add_hline(y=0, line_dash="dash", line_color="gray", line_width=1)
    
    fig.update_layout(
        title=f"Residual Difference (Pol2 - Pol1) vs Max Stimulus Amplitude",
        xaxis_title="Max Stimulus Amplitude (mA)",
        yaxis_title="Residual Difference (mm)",
        hovermode='closest',
        width=1100,
        height=700,
        template='plotly_white',
        font=dict(size=12)
    )
    
    fig.show()


# ── Generate plots ───────────────────────────────────────────────────────────
print("\n" + "="*60)
print("PREDICTIONS USING PEAK AMPLITUDE")
print("="*60)
fit_and_plot_prediction("peak_amplitude", "Peak Amplitude (µV)")
plot_residual_differences("peak_amplitude", "Peak Amplitude (µV)")

print("\n" + "="*60)
print("PREDICTIONS USING PEAK TIME")
print("="*60)
fit_and_plot_prediction("peak_time", "Peak Time (ms)")
plot_residual_differences("peak_time", "Peak Time (ms)")

print("\n=== All plots complete ===")

print("\n=== All plots complete ===")
