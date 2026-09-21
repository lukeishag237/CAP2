import pandas as pd
import plotly.graph_objects as go
import os
import numpy as np

LOCAL_OUTPUT_PATH = "/Users/lukeishag/Documents/PERL/CAP2/output"
df = pd.read_pickle(os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl"))

# Only keep myelinated rows
df = df[df["target_fiber_types"] == "myel"]

# Normalize expdate
df["expdate"] = df["expdate"].astype(str).str.replace(".0", "", regex=False)

# ==================================================================================
# Account for recording cuff contacts before any calculations
# ==================================================================================
def apply_cuff_contact_correction(row):
    """
    Apply recording cuff contact corrections based on nerve type and channel configuration.
    
    For tibial/sciatic:
    - Tripolar (3 channels): subtract 1.5mm from ch1, do nothing to ch2, add 1.5mm to ch3
    - Bipolar (2 channels): subtract 0.75mm from ch1, add 0.75mm to ch2
    
    For cVN (must be tripolar):
    - Tripolar (3 channels): subtract 1.07mm from ch1, do nothing to ch2, add 1.07mm to ch3
    - If only 2 channels: raise error
    """
    nerve = row.get("rec_nerve", "").lower()
    ch_num = row.get("ch_num_neural", 0)
    stim_rec_dist = row.get("stim_rec_dist", 0)
    
    if pd.isna(nerve) or pd.isna(ch_num) or pd.isna(stim_rec_dist):
        return stim_rec_dist
    
    nerve_str = str(nerve).lower()
    
    # Determine if this row represents the maximum channel count
    # We'll check the channel number to infer polarity
    # Note: ch_num_neural is the recording channel
    
    if nerve_str in ["tibial", "sciatic", "sn"]:
        # Check if this appears to be tripolar or bipolar based on the data structure
        # In the dataframe, if there are 3 unique channels, it's tripolar; if 2, it's bipolar
        # For now, we need to infer from ch_num_neural
        # ch_num_neural values: 1, 2, 3 for tripolar; 1, 2 for bipolar
        
        if ch_num == 1:
            return stim_rec_dist - 0.75  # Will be adjusted to 1.5 if tripolar
        elif ch_num == 2:
            return stim_rec_dist + 0.75  # Will be adjusted based on tripolar/bipolar
        elif ch_num == 3:
            # This is tripolar, so go back and fix ch1 and ch2
            return stim_rec_dist + 1.5
        else:
            return stim_rec_dist
    
    elif nerve_str == "cvn":
        if ch_num == 1:
            return stim_rec_dist - 1.07
        elif ch_num == 2:
            return stim_rec_dist  # Do nothing
        elif ch_num == 3:
            return stim_rec_dist + 1.07
        elif ch_num < 3:
            raise ValueError(f"cVN nerve must have 3 channels for proper recording cuff contact correction, but found only {ch_num} channel(s)")
        else:
            return stim_rec_dist
    
    return stim_rec_dist

# First, determine if each nerve type is bipolar or tripolar by checking max channels in the data
def get_nerve_polarity(df):
    """Determine if each nerve type is bipolar or tripolar based on max channels"""
    polarity_map = {}
    for nerve in df["rec_nerve"].unique():
        max_ch = df[df["rec_nerve"] == nerve]["ch_num_neural"].max()
        polarity_map[nerve] = "tripolar" if max_ch >= 3 else "bipolar"
    return polarity_map

polarity_map = get_nerve_polarity(df)

# Apply corrections with proper bipolar/tripolar handling
def apply_cuff_contact_correction_v2(row):
    """Apply recording cuff contact corrections with bipolar/tripolar awareness"""
    nerve = row.get("rec_nerve", "")
    ch_num = row.get("ch_num_neural", 0)
    stim_rec_dist = row.get("stim_rec_dist", 0)
    
    if pd.isna(nerve) or pd.isna(ch_num) or pd.isna(stim_rec_dist):
        return stim_rec_dist
    
    nerve_str = str(nerve).lower()
    is_tripolar = polarity_map.get(nerve, "bipolar") == "tripolar"
    
    if nerve_str in ["tibial", "sciatic", "sn"]:
        if is_tripolar:
            if ch_num == 1:
                return stim_rec_dist - 1.5
            elif ch_num == 2:
                return stim_rec_dist
            elif ch_num == 3:
                return stim_rec_dist + 1.5
        else:  # bipolar
            if ch_num == 1:
                return stim_rec_dist - 0.75
            elif ch_num == 2:
                return stim_rec_dist + 0.75
        return stim_rec_dist
    
    elif nerve_str == "cvn":
        # if ch_num < 3:
        #     print(nerve_str)
        #     raise ValueError(f"cVN nerve must have 3 channels for proper recording cuff contact correction, but found only {ch_num} channel(s)")
        if ch_num == 1:
            return stim_rec_dist - 1.07
        elif ch_num == 2:
            return stim_rec_dist
        elif ch_num == 3:
            return stim_rec_dist + 1.07
    
    return stim_rec_dist

# Apply the cuff contact correction to stim_rec_dist
df["stim_rec_dist"] = df.apply(apply_cuff_contact_correction_v2, axis=1)

# Keys that define a CAP instance
keys = ["expdate", "rec_nerve", "target_fiber_types",
        "ch_num_neural", "stim_rec_dist"]

# Pivot polarity 1 and 2 into columns
pivot = df.pivot_table(
    index=keys,
    columns="dataset_stim_polarity",
    values=["auc", "peak_time", "peak_amplitude", "stim_PW", "max_stim_amplitude_mA"],
    aggfunc="first"
)

# Flatten MultiIndex columns
pivot.columns = [f"{metric}_pol{int(pol)}" for metric, pol in pivot.columns]

# Drop incomplete pairs
pivot = pivot.dropna()

# Build per-pair categorical labels
pivot_reset = pivot.reset_index()

pivot_reset["pair_id"] = (
    pivot_reset["expdate"].astype(str)
    + "_"
    + pivot_reset["rec_nerve"].astype(str)
    + "_ch"
    + pivot_reset["ch_num_neural"].astype(str)
    + "_d"
    + pivot_reset["stim_rec_dist"].astype(str)
)

# Sort by nerve type
nerve_order = {"cVN": 0, "tibial": 1, "SN": 2}
pivot_reset["nerve_sort"] = pivot_reset["rec_nerve"].map(nerve_order)
pivot_reset = pivot_reset.sort_values(["nerve_sort", "expdate"])

# Color map for nerves
nerve_colors = {
    "SN": "#2ca02c",      # Green
    "tibial": "#1f77b4",  # Blue
    "cVN": "#d62728"      # Red
}

# Symbol mapping for channels
CHANNEL_SYMBOLS = {
    1: "circle",
    2: "square",
    3: "diamond"
}

# Calculate actual vs calculated time delays
# Actual time delay = peak_time_pol2 - peak_time_pol1
# (positive if pol2 is delayed, negative if pol1 is delayed)
# Note: polarity 1 = negative leading, polarity 2 = positive leading
# So actual_time_delay = positive_peak_time - negative_peak_time

pivot_reset["actual_time_delay"] = pivot_reset["peak_time_pol2"] - pivot_reset["peak_time_pol1"]
pivot_reset["calculated_time_delay"] = pivot_reset["stim_PW_pol1"]

# Add jitter to help visualize overlapping points
np.random.seed(42)  # For reproducibility
jitter_amount = 0.001  # Small jitter in milliseconds
pivot_reset["actual_time_delay_jitter"] = pivot_reset["actual_time_delay"] + np.random.normal(0, jitter_amount, len(pivot_reset))
pivot_reset["calculated_time_delay_jitter"] = pivot_reset["calculated_time_delay"] + np.random.normal(0, jitter_amount, len(pivot_reset))

# Create Plotly figure
fig = go.Figure()

# Add scatter points for each nerve and channel
for nerve in ["SN", "tibial", "cVN"]:
    nerve_data = pivot_reset[pivot_reset["rec_nerve"] == nerve]
    
    for ch in sorted(nerve_data["ch_num_neural"].unique()):
        ch_data = nerve_data[nerve_data["ch_num_neural"] == ch]
        
        if len(ch_data) == 0:
            continue
        
        symbol = CHANNEL_SYMBOLS.get(ch, "circle")
        
        # Create hover text
        hover_text = [f"Expdate: {row['expdate']}<br>"
                     f"Nerve: {row['rec_nerve']}<br>"
                     f"Channel: {int(row['ch_num_neural'])}<br>"
                     f"Distance: {row['stim_rec_dist']:.2f} mm<br>"
                     f"Actual delay: {row['actual_time_delay']:.3f} ms<br>"
                     f"Calculated delay: {row['calculated_time_delay']:.3f} ms"
                     for _, row in ch_data.iterrows()]
        
        fig.add_trace(go.Scatter(
            x=ch_data["actual_time_delay_jitter"],
            y=ch_data["calculated_time_delay_jitter"],
            mode='markers',
            marker=dict(
                symbol=symbol,
                size=12,
                color=nerve_colors[nerve],
                line=dict(width=1.5, color=nerve_colors[nerve]),
                opacity=0.6  # Reduced opacity to show overlaps
            ),
            name=f"{nerve} Ch{ch}",
            legendgroup=nerve,
            showlegend=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))

# Add diagonal reference line (perfect agreement)
all_x = pivot_reset["actual_time_delay"].dropna()
all_y = pivot_reset["calculated_time_delay"].dropna()
if len(all_x) > 0 and len(all_y) > 0:
    lim_min = min(all_x.min(), all_y.min())
    lim_max = max(all_x.max(), all_y.max())
    padding = (lim_max - lim_min) * 0.05
    lim = [lim_min - padding, lim_max + padding]
    
    fig.add_trace(go.Scatter(
        x=lim,
        y=lim,
        mode='lines',
        line=dict(color='gray', width=2, dash='dash'),
        name='Perfect agreement',
        hoverinfo='skip',
        showlegend=True
    ))

fig.update_layout(
    title="Hypothesis 1 (Front Contact): Actual vs Calculated Time Delay Between Polarities (All Channels, Myelinated Only)",
    xaxis_title="Actual Time Delay (ms)<br>(peak_time_positive - peak_time_negative)",
    yaxis_title="Calculated Time Delay (ms)<br>(stim_PW)",
    hovermode='closest',
    width=1000,
    height=800,
    template='plotly_white',
    font=dict(size=12),
    showlegend=True
)

fig.show()

print(f"Total pairs plotted: {len(pivot_reset)}")
print(f"Pairs per nerve:")
for nerve in ["SN", "tibial", "cVN"]:
    count = len(pivot_reset[pivot_reset["rec_nerve"] == nerve])
    print(f"  {nerve}: {count}")

# ==================================================================================
# HYPOTHESIS 2: Back Contact
# ==================================================================================
# Assumption: Positive leading phase (pol 2) originates 1.5 mm back from negative
# 
# Calculation:
# 1. Calculate CV from negative phase: CV = stim_rec_dist / peak_time_pol1
# 2. Calculate expected peak time for positive phase: 
#    peak_time_expected_pol2 = (stim_rec_dist + 1.5) / CV
# 3. Calculate expected time delay: 
#    calculated_delay_h2 = peak_time_expected_pol2 - peak_time_pol1

print("\n" + "="*80)
print("HYPOTHESIS 2: Back Contact (positive phase originates 1.5mm back)")
print("="*80)

# Calculate conduction velocities and back-solved distances
pivot_reset["CV_ms_per_mm"] = pivot_reset["stim_rec_dist"] / pivot_reset["peak_time_pol1"]
pivot_reset["peak_time_expected_pol2_h2"] = (pivot_reset["stim_rec_dist"] + 1.5) / pivot_reset["CV_ms_per_mm"]
pivot_reset["calculated_time_delay_h2"] = pivot_reset["peak_time_expected_pol2_h2"] - pivot_reset["peak_time_pol1"]

# No jitter for H2 - use raw values

# Create second Plotly figure
fig2 = go.Figure()

# Add scatter points for each nerve and channel
for nerve in ["SN", "tibial", "cVN"]:
    nerve_data = pivot_reset[pivot_reset["rec_nerve"] == nerve]
    
    for ch in sorted(nerve_data["ch_num_neural"].unique()):
        ch_data = nerve_data[nerve_data["ch_num_neural"] == ch]
        
        if len(ch_data) == 0:
            continue
        
        symbol = CHANNEL_SYMBOLS.get(ch, "circle")
        
        # Create hover text
        hover_text = [f"Expdate: {row['expdate']}<br>"
                     f"Nerve: {row['rec_nerve']}<br>"
                     f"Channel: {int(row['ch_num_neural'])}<br>"
                     f"Distance: {row['stim_rec_dist']:.2f} mm<br>"
                     f"CV: {row['CV_ms_per_mm']:.2f} m/s<br>"
                     f"Actual delay: {row['actual_time_delay']:.3f} ms<br>"
                     f"Calculated delay (H2): {row['calculated_time_delay_h2']:.3f} ms"
                     for _, row in ch_data.iterrows()]
        
        fig2.add_trace(go.Scatter(
            x=ch_data["actual_time_delay"],
            y=ch_data["calculated_time_delay_h2"],
            mode='markers',
            marker=dict(
                symbol=symbol,
                size=12,
                color=nerve_colors[nerve],
                line=dict(width=1.5, color=nerve_colors[nerve]),
                opacity=0.6
            ),
            name=f"{nerve} Ch{ch}",
            legendgroup=nerve,
            showlegend=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))

# Add diagonal reference line (perfect agreement)
all_x = pivot_reset["actual_time_delay"].dropna()
all_y = pivot_reset["calculated_time_delay_h2"].dropna()
if len(all_x) > 0 and len(all_y) > 0:
    lim_min = min(all_x.min(), all_y.min())
    lim_max = max(all_x.max(), all_y.max())
    padding = (lim_max - lim_min) * 0.05
    lim = [lim_min - padding, lim_max + padding]
    
    fig2.add_trace(go.Scatter(
        x=lim,
        y=lim,
        mode='lines',
        line=dict(color='gray', width=2, dash='dash'),
        name='Perfect agreement',
        hoverinfo='skip',
        showlegend=True
    ))

fig2.update_layout(
    title="Hypothesis 2 (Back Contact): Actual vs Calculated Time Delay Between Polarities (All Channels, Myelinated Only)",
    xaxis_title="Actual Time Delay (ms)<br>(peak_time_positive - peak_time_negative)",
    yaxis_title="Calculated Time Delay (ms)<br>(back-solved from CV, assuming 1.5mm offset)",
    hovermode='closest',
    width=1000,
    height=800,
    template='plotly_white',
    font=dict(size=12),
    showlegend=True
)

fig2.show()

print(f"Total pairs in H2: {len(pivot_reset)}")
print(f"Pairs per nerve:")
for nerve in ["SN", "tibial", "cVN"]:
    count = len(pivot_reset[pivot_reset["rec_nerve"] == nerve])
    print(f"  {nerve}: {count}")

# ==================================================================================
# HYPOTHESIS 3: Virtual Contact
# ==================================================================================
# X-axis: Maximum stimulus amplitude
# Y-axis: Actual time delay (peak_time_pol2 - peak_time_pol1)

print("\n" + "="*80)
print("HYPOTHESIS 3: Virtual Contact")
print("="*80)

# Get max stim amplitude (should be same for both polarities, but take pol1)
pivot_reset["max_stim_amplitude"] = pivot_reset["max_stim_amplitude_mA_pol1"]

# No jitter for H3 - use raw values

# Create third Plotly figure
fig3 = go.Figure()

# Add scatter points for each nerve and channel
for nerve in ["SN", "tibial", "cVN"]:
    nerve_data = pivot_reset[pivot_reset["rec_nerve"] == nerve]
    
    for ch in sorted(nerve_data["ch_num_neural"].unique()):
        ch_data = nerve_data[nerve_data["ch_num_neural"] == ch]
        
        if len(ch_data) == 0:
            continue
        
        symbol = CHANNEL_SYMBOLS.get(ch, "circle")
        
        # Create hover text
        hover_text = [f"Expdate: {row['expdate']}<br>"
                     f"Nerve: {row['rec_nerve']}<br>"
                     f"Channel: {int(row['ch_num_neural'])}<br>"
                     f"Distance: {row['stim_rec_dist']:.2f} mm<br>"
                     f"Max Stimulus Amplitude: {row['max_stim_amplitude']:.3f} mA<br>"
                     f"Actual time delay: {row['actual_time_delay']:.3f} ms"
                     for _, row in ch_data.iterrows()]
        
        fig3.add_trace(go.Scatter(
            x=ch_data["max_stim_amplitude"],
            y=ch_data["actual_time_delay"],
            mode='markers',
            marker=dict(
                symbol=symbol,
                size=12,
                color=nerve_colors[nerve],
                line=dict(width=1.5, color=nerve_colors[nerve]),
                opacity=0.6
            ),
            name=f"{nerve} Ch{ch}",
            legendgroup=nerve,
            showlegend=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))

fig3.update_layout(
    title="Hypothesis 3 (Virtual Contact): Actual Time Delay vs Stimulus Amplitude (All Channels, Myelinated Only)",
    xaxis_title="Maximum Stimulus Amplitude (mA)",
    yaxis_title="Actual Time Delay (ms)<br>(peak_time_positive - peak_time_negative)",
    hovermode='closest',
    width=1000,
    height=800,
    template='plotly_white',
    font=dict(size=12),
    showlegend=True
)

fig3.show()

print(f"Total pairs in H3: {len(pivot_reset)}")
print(f"Pairs per nerve:")
for nerve in ["SN", "tibial", "cVN"]:
    count = len(pivot_reset[pivot_reset["rec_nerve"] == nerve])
    print(f"  {nerve}: {count}")

# ==================================================================================
# HYPOTHESIS 4: Combined (H1 + H2)
# ==================================================================================
# Calculate combined calculated delay: (stim_PW + back_solved_delay) / 2
# Or add them together: calculated_delay_h1 + calculated_delay_h2

print("\n" + "="*80)
print("HYPOTHESIS 4: Combined (H1 + H2)")
print("="*80)

# Combine the calculated delays from H1 and H2
pivot_reset["calculated_time_delay_h4"] = pivot_reset["calculated_time_delay"] + pivot_reset["calculated_time_delay_h2"]

# Create fourth Plotly figure
fig4 = go.Figure()

# Add scatter points for each nerve and channel
for nerve in ["SN", "tibial", "cVN"]:
    nerve_data = pivot_reset[pivot_reset["rec_nerve"] == nerve]
    
    for ch in sorted(nerve_data["ch_num_neural"].unique()):
        ch_data = nerve_data[nerve_data["ch_num_neural"] == ch]
        
        if len(ch_data) == 0:
            continue
        
        symbol = CHANNEL_SYMBOLS.get(ch, "circle")
        
        # Create hover text
        hover_text = [f"Expdate: {row['expdate']}<br>"
                     f"Nerve: {row['rec_nerve']}<br>"
                     f"Channel: {int(row['ch_num_neural'])}<br>"
                     f"Distance: {row['stim_rec_dist']:.2f} mm<br>"
                     f"Actual delay: {row['actual_time_delay']:.3f} ms<br>"
                     f"H1 (stim_PW): {row['calculated_time_delay']:.3f} ms<br>"
                     f"H2 (back-solved): {row['calculated_time_delay_h2']:.3f} ms<br>"
                     f"H4 (H1+H2): {row['calculated_time_delay_h4']:.3f} ms"
                     for _, row in ch_data.iterrows()]
        
        fig4.add_trace(go.Scatter(
            x=ch_data["actual_time_delay"],
            y=ch_data["calculated_time_delay_h4"],
            mode='markers',
            marker=dict(
                symbol=symbol,
                size=12,
                color=nerve_colors[nerve],
                line=dict(width=1.5, color=nerve_colors[nerve]),
                opacity=0.6
            ),
            name=f"{nerve} Ch{ch}",
            legendgroup=nerve,
            showlegend=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))

# Add diagonal reference line (perfect agreement)
all_x = pivot_reset["actual_time_delay"].dropna()
all_y = pivot_reset["calculated_time_delay_h4"].dropna()
if len(all_x) > 0 and len(all_y) > 0:
    lim_min = min(all_x.min(), all_y.min())
    lim_max = max(all_x.max(), all_y.max())
    padding = (lim_max - lim_min) * 0.05
    lim = [lim_min - padding, lim_max + padding]
    
    fig4.add_trace(go.Scatter(
        x=lim,
        y=lim,
        mode='lines',
        line=dict(color='gray', width=2, dash='dash'),
        name='Perfect agreement',
        hoverinfo='skip',
        showlegend=True
    ))

fig4.update_layout(
    title="Hypothesis 4 (Combined H1+H2): Actual vs Combined Calculated Time Delay (All Channels, Myelinated Only)",
    xaxis_title="Actual Time Delay (ms)<br>(peak_time_positive - peak_time_negative)",
    yaxis_title="Combined Calculated Time Delay (ms)<br>(stim_PW + back-solved)",
    hovermode='closest',
    width=1000,
    height=800,
    template='plotly_white',
    font=dict(size=12),
    showlegend=True
)

fig4.show()

print(f"Total pairs in H4: {len(pivot_reset)}")
print(f"Pairs per nerve:")
for nerve in ["SN", "tibial", "cVN"]:
    count = len(pivot_reset[pivot_reset["rec_nerve"] == nerve])
    print(f"  {nerve}: {count}")
