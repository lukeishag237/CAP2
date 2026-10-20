import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os
import numpy as np

# Load data
script_dir = os.path.dirname(os.path.abspath(__file__))
actual_df = pd.read_csv(os.path.join(script_dir, "actual_values.csv"))
paired_df = pd.read_csv(os.path.join(script_dir, "paired_difference.csv"))

# Nerve type colors
NERVE_COLORS = {
    "sciatic": "#d62728",  # Red
    "tibial": "#1f77b4"    # Blue
}

# Polarity shapes
POLARITY_SYMBOLS = {
    1: "circle",
    2: "square"
}

print("Data loaded successfully!")
print(f"Actual values shape: {actual_df.shape}")
print(f"Paired differences shape: {paired_df.shape}")
print(f"\nActual values columns: {list(actual_df.columns)}")
print(f"Paired differences columns: {list(paired_df.columns)}")

# ==================================================================================
# PREPARE DATA: Expand by channel and apply CD corrections
# ==================================================================================

def expand_and_correct_cd(df, is_paired=False):
    """
    Expand dataframe so each channel is a separate row, and apply CD corrections.
    
    For tripolar (3 channels): subtract 1.5mm from ch1, leave ch2, add 1.5mm to ch3
    For bipolar (2 channels): subtract 0.75mm from ch1, add 0.75mm to ch2
    """
    expanded_rows = []
    
    for idx, row in df.iterrows():
        # Check how many channels we have by looking for NA values
        has_ch1 = not pd.isna(row.get("Ch1 time"))
        has_ch2 = not pd.isna(row.get("Ch2 time"))
        has_ch3 = not pd.isna(row.get("Ch3 time"))
        
        # Determine if tripolar or bipolar
        is_tripolar = has_ch3
        
        # Create a row for each channel
        if has_ch1:
            ch1_row = row.copy()
            ch1_row["Channel"] = 1
            if is_tripolar:
                ch1_row["Adjusted_CD"] = row["Conduction Distance (mm)"] - 1.5
            else:
                ch1_row["Adjusted_CD"] = row["Conduction Distance (mm)"] - 0.75
            ch1_row["Time"] = row.get("Ch1 time")
            ch1_row["P2P"] = row.get("Ch1 p2p")
            ch1_row["Spread"] = row.get("Ch1 spread")
            expanded_rows.append(ch1_row)
        
        if has_ch2:
            ch2_row = row.copy()
            ch2_row["Channel"] = 2
            if is_tripolar:
                ch2_row["Adjusted_CD"] = row["Conduction Distance (mm)"]  # No change
            else:
                ch2_row["Adjusted_CD"] = row["Conduction Distance (mm)"] + 0.75
            ch2_row["Time"] = row.get("Ch2 time")
            ch2_row["P2P"] = row.get("Ch2 p2p")
            ch2_row["Spread"] = row.get("Ch2 spread")
            expanded_rows.append(ch2_row)
        
        if has_ch3:
            ch3_row = row.copy()
            ch3_row["Channel"] = 3
            ch3_row["Adjusted_CD"] = row["Conduction Distance (mm)"] + 1.5
            ch3_row["Time"] = row.get("Ch3 time")
            ch3_row["P2P"] = row.get("Ch3 p2p")
            ch3_row["Spread"] = row.get("Ch3 spread")
            expanded_rows.append(ch3_row)
    
    expanded_df = pd.DataFrame(expanded_rows)
    return expanded_df

# Expand both dataframes
actual_expanded = expand_and_correct_cd(actual_df, is_paired=False)
paired_expanded = expand_and_correct_cd(paired_df, is_paired=True)

print(f"\nExpanded actual values shape: {actual_expanded.shape}")
print(f"Expanded paired differences shape: {paired_expanded.shape}")
print(f"Nerve types: {actual_expanded['Nerve type'].unique()}")

# ==================================================================================
# ANALYSIS 1: TIME MEASUREMENTS (vs Adjusted CD, PW, and Stim Amplitude)
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 1: TIME MEASUREMENTS")
print("="*80)

fig1 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("Time vs Adjusted Conduction Distance", "Time vs Pulse Width", "Time vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# Time vs Adjusted Conduction Distance
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                     f"Original CD: {row['Conduction Distance (mm)']:.2f} mm<br>"
                     f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                     f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                     f"Time: {row['Time']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig1.add_trace(
            go.Scatter(
                x=nerve_pol_data["Adjusted_CD"],
                y=nerve_pol_data["Time"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=True
            ),
            row=1, col=1
        )

# Time vs Pulse Width
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                     f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                     f"Time: {row['Time']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig1.add_trace(
            go.Scatter(
                x=nerve_pol_data["Pulse Width (ms)"],
                y=nerve_pol_data["Time"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=2
        )

# Time vs Stim Amplitude
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                     f"Time: {row['Time']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig1.add_trace(
            go.Scatter(
                x=nerve_pol_data["Stim Amplitude (mA)"],
                y=nerve_pol_data["Time"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=3
        )

fig1.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig1.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig1.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig1.update_yaxes(title_text="Time (ms)", row=1, col=1)

fig1.update_layout(
    title="Actual Values: Time Measurements vs Parameters (Each Channel as Separate Measurement)",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig1.show()

# ==================================================================================
# ANALYSIS 2: PEAK-TO-PEAK AMPLITUDE
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 2: PEAK-TO-PEAK AMPLITUDE")
print("="*80)

fig2 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("P2P vs Adjusted Conduction Distance", "P2P vs Pulse Width", "P2P vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# P2P vs Adjusted Conduction Distance
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                     f"P2P: {row['P2P']:.2f} µV<br>"
                     f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig2.add_trace(
            go.Scatter(
                x=nerve_pol_data["Adjusted_CD"],
                y=nerve_pol_data["P2P"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=1
        )

# P2P vs Pulse Width
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                     f"P2P: {row['P2P']:.2f} µV"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig2.add_trace(
            go.Scatter(
                x=nerve_pol_data["Pulse Width (ms)"],
                y=nerve_pol_data["P2P"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=2
        )

# P2P vs Stim Amplitude
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                     f"P2P: {row['P2P']:.2f} µV"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig2.add_trace(
            go.Scatter(
                x=nerve_pol_data["Stim Amplitude (mA)"],
                y=nerve_pol_data["P2P"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=3
        )

fig2.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig2.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig2.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig2.update_yaxes(title_text="Peak-to-Peak Amplitude (µV)", row=1, col=1)

fig2.update_layout(
    title="Actual Values: Peak-to-Peak Amplitude vs Parameters",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig2.show()

# ==================================================================================
# ANALYSIS 3: SIGNAL SPREAD
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 3: SIGNAL SPREAD")
print("="*80)

fig3 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("Spread vs Adjusted Conduction Distance", "Spread vs Pulse Width", "Spread vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# Spread vs Adjusted Conduction Distance
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                     f"Spread: {row['Spread']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig3.add_trace(
            go.Scatter(
                x=nerve_pol_data["Adjusted_CD"],
                y=nerve_pol_data["Spread"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=1
        )

# Spread vs Pulse Width
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                     f"Spread: {row['Spread']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig3.add_trace(
            go.Scatter(
                x=nerve_pol_data["Pulse Width (ms)"],
                y=nerve_pol_data["Spread"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=2
        )

# Spread vs Stim Amplitude
for nerve in sorted(actual_expanded["Nerve type"].unique()):
    for polarity in sorted(actual_expanded["Polarity"].unique()):
        nerve_pol_data = actual_expanded[(actual_expanded["Nerve type"] == nerve) & 
                                         (actual_expanded["Polarity"] == polarity)]
        if len(nerve_pol_data) == 0:
            continue
        
        hover_text = [f"ID: {row['Exp ID']}<br>"
                     f"Nerve: {row['Nerve type']}<br>"
                     f"Channel: {int(row['Channel'])}<br>"
                     f"Polarity: {int(row['Polarity'])}<br>"
                     f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                     f"Spread: {row['Spread']:.3f} ms"
                     for _, row in nerve_pol_data.iterrows()]
        
        fig3.add_trace(
            go.Scatter(
                x=nerve_pol_data["Stim Amplitude (mA)"],
                y=nerve_pol_data["Spread"],
                mode='markers',
                name=f"{nerve} Pol{int(polarity)}",
                marker=dict(
                    color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                    symbol=POLARITY_SYMBOLS.get(polarity, "circle"),
                    size=10,
                    opacity=0.7,
                    line=dict(width=1)
                ),
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False,
                legendgroup=f"{nerve}_{polarity}"
            ),
            row=1, col=3
        )

fig3.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig3.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig3.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig3.update_yaxes(title_text="Signal Spread (ms)", row=1, col=1)

fig3.update_layout(
    title="Actual Values: Signal Spread vs Parameters",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig3.show()

# ==================================================================================
# ANALYSIS 4: PAIRED DIFFERENCES - TIME
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 4: PAIRED DIFFERENCES - TIME")
print("="*80)

fig4 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("Time Diff vs Adjusted Conduction Distance", "Time Diff vs Pulse Width", "Time Diff vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# Time differences vs Adjusted Conduction Distance
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                 f"Time Diff (Pol2-Pol1): {row['Time']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig4.add_trace(
        go.Scatter(
            x=nerve_data["Adjusted_CD"],
            y=nerve_data["Time"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True
        ),
        row=1, col=1
    )

# Time differences vs Pulse Width
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                 f"Time Diff (Pol2-Pol1): {row['Time']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig4.add_trace(
        go.Scatter(
            x=nerve_data["Pulse Width (ms)"],
            y=nerve_data["Time"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=2
    )

# Time differences vs Stim Amplitude
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                 f"Time Diff (Pol2-Pol1): {row['Time']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig4.add_trace(
        go.Scatter(
            x=nerve_data["Stim Amplitude (mA)"],
            y=nerve_data["Time"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=3
    )

# Add horizontal line at y=0 to show no difference
fig4.add_hline(y=0, line_dash="dash", line_color="gray", annotation_text="No difference", row="all", col="all")

fig4.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig4.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig4.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig4.update_yaxes(title_text="Time Difference (ms)", row=1, col=1)

fig4.update_layout(
    title="Paired Differences: Time Differences vs Parameters (Pol2 - Pol1)",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig4.show()

# ==================================================================================
# ANALYSIS 5: PAIRED DIFFERENCES - P2P
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 5: PAIRED DIFFERENCES - P2P")
print("="*80)

fig5 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("P2P Diff vs Adjusted Conduction Distance", "P2P Diff vs Pulse Width", "P2P Diff vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# P2P differences vs Adjusted Conduction Distance
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                 f"P2P Diff (Pol2-Pol1): {row['P2P']:.2f} µV"
                 for _, row in nerve_data.iterrows()]
    
    fig5.add_trace(
        go.Scatter(
            x=nerve_data["Adjusted_CD"],
            y=nerve_data["P2P"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True
        ),
        row=1, col=1
    )

# P2P differences vs Pulse Width
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                 f"P2P Diff (Pol2-Pol1): {row['P2P']:.2f} µV"
                 for _, row in nerve_data.iterrows()]
    
    fig5.add_trace(
        go.Scatter(
            x=nerve_data["Pulse Width (ms)"],
            y=nerve_data["P2P"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=2
    )

# P2P differences vs Stim Amplitude
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                 f"P2P Diff (Pol2-Pol1): {row['P2P']:.2f} µV"
                 for _, row in nerve_data.iterrows()]
    
    fig5.add_trace(
        go.Scatter(
            x=nerve_data["Stim Amplitude (mA)"],
            y=nerve_data["P2P"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=3
    )

# Add horizontal line at y=0 to show no difference
fig5.add_hline(y=0, line_dash="dash", line_color="gray", annotation_text="No difference", row="all", col="all")

fig5.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig5.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig5.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig5.update_yaxes(title_text="P2P Difference (µV)", row=1, col=1)

fig5.update_layout(
    title="Paired Differences: Peak-to-Peak Amplitude Differences vs Parameters (Pol2 - Pol1)",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig5.show()

# ==================================================================================
# ANALYSIS 6: PAIRED DIFFERENCES - SPREAD
# ==================================================================================
print("\n" + "="*80)
print("ANALYSIS 6: PAIRED DIFFERENCES - SPREAD")
print("="*80)

fig6 = make_subplots(
    rows=1, cols=3,
    subplot_titles=("Spread Diff vs Adjusted Conduction Distance", "Spread Diff vs Pulse Width", "Spread Diff vs Stim Amplitude"),
    specs=[[{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}]]
)

# Spread differences vs Adjusted Conduction Distance
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Adjusted CD: {row['Adjusted_CD']:.2f} mm<br>"
                 f"Spread Diff (Pol2-Pol1): {row['Spread']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig6.add_trace(
        go.Scatter(
            x=nerve_data["Adjusted_CD"],
            y=nerve_data["Spread"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=True
        ),
        row=1, col=1
    )

# Spread differences vs Pulse Width
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Pulse Width: {row['Pulse Width (ms)']:.4f} ms<br>"
                 f"Spread Diff (Pol2-Pol1): {row['Spread']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig6.add_trace(
        go.Scatter(
            x=nerve_data["Pulse Width (ms)"],
            y=nerve_data["Spread"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=2
    )

# Spread differences vs Stim Amplitude
for nerve in sorted(paired_expanded["Nerve type"].unique()):
    nerve_data = paired_expanded[paired_expanded["Nerve type"] == nerve]
    if len(nerve_data) == 0:
        continue
    
    hover_text = [f"Pair ID: {row['Pair IDs']}<br>"
                 f"Nerve: {row['Nerve type']}<br>"
                 f"Channel: {int(row['Channel'])}<br>"
                 f"Stim Amplitude: {row['Stim Amplitude (mA)']:.2f} mA<br>"
                 f"Spread Diff (Pol2-Pol1): {row['Spread']:.3f} ms"
                 for _, row in nerve_data.iterrows()]
    
    fig6.add_trace(
        go.Scatter(
            x=nerve_data["Stim Amplitude (mA)"],
            y=nerve_data["Spread"],
            mode='markers',
            name=nerve,
            marker=dict(
                color=NERVE_COLORS.get(nerve.lower(), "#000000"),
                symbol="circle",
                size=10,
                opacity=0.7,
                line=dict(width=1)
            ),
            hovertext=hover_text,
            hoverinfo='text',
            showlegend=False
        ),
        row=1, col=3
    )

# Add horizontal line at y=0 to show no difference
fig6.add_hline(y=0, line_dash="dash", line_color="gray", annotation_text="No difference", row="all", col="all")

fig6.update_xaxes(title_text="Adjusted Conduction Distance (mm)", row=1, col=1)
fig6.update_xaxes(title_text="Pulse Width (ms)", row=1, col=2)
fig6.update_xaxes(title_text="Stim Amplitude (mA)", row=1, col=3)
fig6.update_yaxes(title_text="Spread Difference (ms)", row=1, col=1)

fig6.update_layout(
    title="Paired Differences: Signal Spread Differences vs Parameters (Pol2 - Pol1)",
    height=500,
    width=1400,
    template="plotly_white",
    hovermode="closest",
    showlegend=True
)

fig6.show()

print("\n" + "="*80)
print("Analysis complete!")
print("="*80)
