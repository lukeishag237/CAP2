#!/usr/bin/env python3
"""
plot_polarity_comparison.py

Uses stored_df.pkl (all myelinated recordings) to produce 4 figures:

  Fig 1 — Unpaired: Peak-to-Peak vs adjusted CD
           colour = nerve (SN/tibial/cVN), shape = channel (1/2/3)
           arrows pol1 → pol2 per (expdate, dist, stim_PW, stim_amp, ch)

  Fig 2 — Unpaired: Peak Time vs adjusted CD  (same encoding)

  Fig 3 — Paired differences, Peak Time (4-panel):
           col1: Pol1 − Pol2  |  col2: Ch1 − Ch2  |  col3: Ch1 − Ch3  |  col4: Ch2 − Ch3

  Fig 4 — Paired differences, Peak-to-Peak (4-panel, same structure)
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ──────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ──────────────────────────────────────────────────────────────────────────────
STORED_DF_PATH = '/Users/lukeishag/Documents/PERL/CAP2/output/stored_df.pkl'
EXCLUDED_PATH  = '/Users/lukeishag/Documents/PERL/CAP2/excluded_datasets.csv'
EXCEL_PATH     = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/raw_data/AnalysisScript_CAPs_Batch_parameters.xlsx'

# Visual constants
NERVE_COLORS = {
    'sn':     '#d62728',   # Red   — Sciatic
    'tibial': '#1f77b4',   # Blue  — Tibial
    'cvn':    '#2ca02c',   # Green — cVN
}
CHANNEL_SYMBOLS = {1: 'circle', 2: 'square', 3: 'diamond'}

def nerve_color(nerve: str) -> str:
    return NERVE_COLORS.get(nerve.lower().strip(), '#7f7f7f')

def nerve_label(nerve: str) -> str:
    mapping = {'sn': 'Sciatic', 'tibial': 'Tibial', 'cvn': 'cVN'}
    return mapping.get(nerve.lower().strip(), nerve)

# ──────────────────────────────────────────────────────────────────────────────
# LOAD & FILTER DATA
# ──────────────────────────────────────────────────────────────────────────────
print("Loading data...")
raw_df = pd.read_pickle(STORED_DF_PATH)
df = raw_df[raw_df['target_fiber_types'] == 'myel'].copy().reset_index(drop=True)
print(f"  stored_df: {len(raw_df)} total rows → {len(df)} myelinated")

# ── Load excluded list (always) ───────────────────────────────────────────────
excl_set = set()
try:
    excl_df = pd.read_csv(EXCLUDED_PATH)
    excl_set = set(excl_df['trial_unique_ID'].dropna().astype(str).str.strip())
    print(f"  {len(excl_set)} excluded IDs loaded")
except Exception as e:
    print(f"  Warning: could not load excluded list: {e}")

# ── Load Excel (always, separately) ──────────────────────────────────────────
myel_xl = None
try:
    print("  Reading Excel metadata (this may take a moment)...")
    xlsx = pd.read_excel(EXCEL_PATH)
    myel_xl = xlsx[xlsx['target_fiber_types'] == 'myel'].copy()
    print(f"  Excel loaded: {len(myel_xl)} myelinated rows")

    # Build exclusion conditions and filter df
    excl_conditions = set()
    for _, row in myel_xl.iterrows():
        prefix = str(row['file_prefix_outputs_values']).strip()
        if prefix in excl_set:
            excl_conditions.add((
                str(int(row['expdate_values'])),
                float(row['stim_rec_distance_center_to_center_mm']),
                float(row['dataset_stim_polarity']),
            ))

    before = len(df)
    df = df[~df.apply(
        lambda r: (str(r['expdate']), float(r['stim_rec_dist']), float(r['dataset_stim_polarity']))
                  in excl_conditions,
        axis=1
    )].reset_index(drop=True)
    print(f"  Excluded {before - len(df)} rows → {len(df)} remain")

except Exception as e:
    print(f"  Warning: Excel load failed ({e}). Exclusion filtering and missing-data report skipped.")

print(f"  Nerves: {sorted(df['rec_nerve'].unique())}")
print(f"  Polarities: {sorted(df['dataset_stim_polarity'].unique())}")
print(f"  Channels: {sorted(df['ch_num_neural'].unique())}")

# ──────────────────────────────────────────────────────────────────────────────
# CD CORRECTION
# ──────────────────────────────────────────────────────────────────────────────
# Determine tripolar vs bipolar per (expdate, stim_rec_dist, pol, PW, amp):
# if ch3 exists in that condition → tripolar (offsets ±1.5 mm)
# otherwise bipolar (offsets ±0.75 mm)

def compute_adjusted_cd(df: pd.DataFrame) -> pd.Series:
    """Add Adjusted_CD column based on tripolar/bipolar detection per condition."""
    # For each row, check if ch3 exists in the same (expdate, dist, pol, PW, amp) group
    group_keys = ['expdate', 'stim_rec_dist', 'dataset_stim_polarity', 'stim_PW', 'max_stim_amplitude_mA', 'rec_nerve']
    has_ch3 = df.groupby(group_keys)['ch_num_neural'].transform(lambda x: 3 in x.values)

    offsets_tripolar  = {1: -1.5, 2:  0.0, 3: +1.5}
    offsets_bipolar   = {1: -0.75, 2: +0.75}

    adj = df.apply(
        lambda r: r['stim_rec_dist'] + (
            offsets_tripolar.get(r['ch_num_neural'], 0.0) if has_ch3[r.name]
            else offsets_bipolar.get(r['ch_num_neural'], 0.0)
        ),
        axis=1
    )
    return adj

df['Adjusted_CD'] = compute_adjusted_cd(df)
df['is_tripolar'] = df.groupby(
    ['expdate','stim_rec_dist','dataset_stim_polarity','stim_PW','max_stim_amplitude_mA','rec_nerve']
)['ch_num_neural'].transform(lambda x: 3 in x.values)

print("  Adjusted_CD computed.\n")

# ──────────────────────────────────────────────────────────────────────────────
# REPORT: Excel rows not in current plot
# ──────────────────────────────────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────────────────────
# REPORT: Excel rows not in current plot
# ──────────────────────────────────────────────────────────────────────────────
if myel_xl is not None:
    stored_keys = set(
        zip(df['expdate'].astype(str),
            df['stim_rec_dist'].round(4),
            df['dataset_stim_polarity'],
            df['stim_PW'].round(4))
    )
    not_plotted = []
    for _, row in myel_xl.iterrows():
        prefix  = str(row['file_prefix_outputs_values']).strip()
        exp     = str(int(row['expdate_values']))
        dist    = round(float(row['stim_rec_distance_center_to_center_mm']), 4)
        pol     = row['dataset_stim_polarity']
        pw      = round(float(row['stim_PW']), 4)
        is_excl = prefix in excl_set
        in_plot = (exp, dist, pol, pw) in stored_keys if not pd.isna(pol) else False
        if not in_plot:
            reason = 'excluded' if is_excl else ('polarity=NaN' if pd.isna(pol) else 'missing from stored_df')
            not_plotted.append((prefix, exp, row['rec_nerve'], dist, pol, pw, reason))

    print("=" * 70)
    print(f"EXCEL ROWS NOT BEING PLOTTED ({len(not_plotted)} total):")
    print("=" * 70)
    for prefix, exp, nerve, dist, pol, pw, reason in not_plotted:
        tag = '[EXCLUDED]' if reason == 'excluded' else '[MISSING] '
        print(f"  {tag}  {prefix:20s}  {nerve:6s}  dist={dist:5.1f}  pol={str(pol):3s}  PW={pw}  → {reason}")
    print("=" * 70 + "\n")
else:
    print("  Skipping missing-data report (Excel not available).\n")

# ──────────────────────────────────────────────────────────────────────────────
# HELPER: legend-dedup tracker
# ──────────────────────────────────────────────────────────────────────────────
def make_legend_tracker():
    seen = set()
    def should_show(key):
        if key in seen:
            return False
        seen.add(key)
        return True
    return should_show

# ──────────────────────────────────────────────────────────────────────────────
# FIG 1 & 2 — UNPAIRED SCATTER WITH ARROWS
# ──────────────────────────────────────────────────────────────────────────────
def make_unpaired_figure(metric_col: str, y_label: str) -> go.Figure:
    """
    One point per row, coloured by nerve, shaped by channel.
    Arrows from pol1 → pol2 for each matched (expdate, dist, PW, amp, ch) group.
    """
    fig = go.Figure()

    # ── Build pol1/pol2 lookup for arrows ────────────────────────────────
    # Key: (expdate, stim_rec_dist, stim_PW, max_stim_amp, rec_nerve, ch_num_neural)
    match_keys = ['expdate', 'stim_rec_dist', 'stim_PW', 'max_stim_amplitude_mA', 'rec_nerve', 'ch_num_neural']
    pol1 = df[df['dataset_stim_polarity'] == 1.0].copy()
    pol2 = df[df['dataset_stim_polarity'] == 2.0].copy()

    # Aggregate duplicates (repeated trials) by mean for arrow endpoints
    pol1_mean = pol1.groupby(match_keys).agg(
        adj_cd=('Adjusted_CD', 'mean'),
        metric=(metric_col, 'mean')
    ).reset_index()
    pol2_mean = pol2.groupby(match_keys).agg(
        adj_cd=('Adjusted_CD', 'mean'),
        metric=(metric_col, 'mean')
    ).reset_index()

    merged = pol1_mean.merge(pol2_mean, on=match_keys, suffixes=('_p1', '_p2'))

    # ── Draw arrows (pol1 → pol2) ─────────────────────────────────────────
    for _, row in merged.iterrows():
        color = nerve_color(row['rec_nerve'])
        fig.add_annotation(
            x=row['adj_cd_p2'], y=row['metric_p2'],
            ax=row['adj_cd_p1'], ay=row['metric_p1'],
            xref='x', yref='y', axref='x', ayref='y',
            showarrow=True,
            arrowhead=3, arrowsize=1.2, arrowwidth=1.8,
            arrowcolor=color,
            opacity=0.45,
        )

    # ── Plot individual points ─────────────────────────────────────────────
    # Legend: one entry per nerve (colour only).
    # Channel shape and polarity are encoded visually but not repeated in legend.
    nerve_legend_shown = set()

    for nerve in sorted(df['rec_nerve'].unique()):
        for ch in [1, 2, 3]:
            for pol in [1, 2]:
                sub = df[
                    (df['rec_nerve'] == nerve) &
                    (df['ch_num_neural'] == ch) &
                    (df['dataset_stim_polarity'] == pol)
                ]
                if sub.empty:
                    continue

                first_for_nerve = nerve.lower() not in nerve_legend_shown
                if first_for_nerve:
                    nerve_legend_shown.add(nerve.lower())

                hover = [
                    f"Exp: {r['expdate']}<br>"
                    f"Nerve: {nerve_label(nerve)}<br>"
                    f"Ch: {ch} | Pol: {int(pol)}<br>"
                    f"Adj CD: {r['Adjusted_CD']:.2f} mm<br>"
                    f"Orig CD: {r['stim_rec_dist']:.1f} mm<br>"
                    f"PW: {r['stim_PW']:.3f} ms | Amp: {r['max_stim_amplitude_mA']:.2f} mA<br>"
                    f"{y_label}: {r[metric_col]:.4f}"
                    for _, r in sub.iterrows()
                ]

                fig.add_trace(go.Scatter(
                    x=sub['Adjusted_CD'],
                    y=sub[metric_col],
                    mode='markers',
                    name=nerve_label(nerve),
                    legendgroup=nerve.lower(),
                    showlegend=first_for_nerve,
                    marker=dict(
                        color=nerve_color(nerve),
                        symbol=CHANNEL_SYMBOLS.get(ch, 'circle'),
                        size=10,
                        opacity=0.85,
                        line=dict(color='black', width=1),
                    ),
                    hovertext=hover,
                    hoverinfo='text',
                ))

    fig.add_hline(y=0, line_dash='dot', line_color='lightgray', line_width=1)

    fig.update_layout(
        title=f'Unpaired — {y_label} vs Adjusted Conduction Distance<br>'
              f'<sup>Colour = nerve type | Shape = channel | Arrows: Pol1 → Pol2</sup>',
        xaxis_title='Adjusted Conduction Distance (mm)',
        yaxis_title=y_label,
        template='plotly_white',
        hovermode='closest',
        width=950, height=620,
        font=dict(size=13),
        legend=dict(groupclick='toggleitem'),
    )
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# FIG 3 & 4 — PAIRED DIFFERENCE (4-panel)
# ──────────────────────────────────────────────────────────────────────────────
def make_paired_diff_figure(metric_col: str, metric_label: str, y_unit: str) -> go.Figure:
    """
    4 sub-plots vs Adjusted CD:
      col1  Pol1 − Pol2  (matched by expdate, dist, PW, amp, ch)
      col2  Ch1  − Ch2   (per polarity, same expdate/dist/PW/amp/nerve)
      col3  Ch1  − Ch3   (tripolar only)
      col4  Ch2  − Ch3   (tripolar only)
    """
    subtitles = [
        f'{metric_label}: Pol1 − Pol2',
        f'{metric_label}: Ch1 − Ch2',
        f'{metric_label}: Ch1 − Ch3',
        f'{metric_label}: Ch2 − Ch3',
    ]
    fig = make_subplots(rows=1, cols=4, subplot_titles=subtitles)
    # One legend entry per nerve across all panels
    nerve_legend_shown = set()

    match_stim = ['expdate', 'stim_rec_dist', 'stim_PW', 'max_stim_amplitude_mA', 'rec_nerve']

    # ── Col 1: Pol1 − Pol2 ────────────────────────────────────────────────
    pol1 = df[df['dataset_stim_polarity'] == 1.0].copy()
    pol2 = df[df['dataset_stim_polarity'] == 2.0].copy()

    merged_pol = pol1.merge(
        pol2,
        on=match_stim + ['ch_num_neural'],
        suffixes=('_p1', '_p2')
    )
    merged_pol['Diff'] = merged_pol[f'{metric_col}_p1'] - merged_pol[f'{metric_col}_p2']
    merged_pol['Adj_CD'] = (merged_pol['Adjusted_CD_p1'] + merged_pol['Adjusted_CD_p2']) / 2.0

    for nerve in sorted(merged_pol['rec_nerve'].unique()):
        sub = merged_pol[merged_pol['rec_nerve'] == nerve]
        first_for_nerve = nerve.lower() not in nerve_legend_shown
        if first_for_nerve:
            nerve_legend_shown.add(nerve.lower())

        hover = [
            f"Exp: {r['expdate']}<br>"
            f"Nerve: {nerve_label(nerve)}<br>"
            f"Ch: {int(r['ch_num_neural'])}<br>"
            f"Adj CD: {r['Adj_CD']:.2f} mm<br>"
            f"PW: {r['stim_PW']:.3f} ms | Amp: {r['max_stim_amplitude_mA']:.2f} mA<br>"
            f"Pol1={r[f'{metric_col}_p1']:.4f}, Pol2={r[f'{metric_col}_p2']:.4f}<br>"
            f"Diff: {r['Diff']:.4f} {y_unit}"
            for _, r in sub.iterrows()
        ]

        fig.add_trace(go.Scatter(
            x=sub['Adj_CD'], y=sub['Diff'],
            mode='markers',
            name=nerve_label(nerve),
            legendgroup=nerve.lower(),
            showlegend=first_for_nerve,
            marker=dict(
                color=nerve_color(nerve), symbol='circle',
                size=9, opacity=0.8,
                line=dict(color='black', width=1),
            ),
            hovertext=hover, hoverinfo='text',
        ), row=1, col=1)

    # ── Cols 2–4: channel differences ─────────────────────────────────────
    ch_pairs = [(1, 2), (1, 3), (2, 3)]
    match_cond = match_stim + ['dataset_stim_polarity']

    for col_idx, (ch_a, ch_b) in enumerate(ch_pairs, start=2):
        sub_a = df[df['ch_num_neural'] == ch_a].copy()
        sub_b = df[df['ch_num_neural'] == ch_b].copy()

        merged_ch = sub_a.merge(
            sub_b,
            on=match_cond,
            suffixes=(f'_ch{ch_a}', f'_ch{ch_b}')
        )
        if merged_ch.empty:
            continue

        merged_ch['Diff'] = merged_ch[f'{metric_col}_ch{ch_a}'] - merged_ch[f'{metric_col}_ch{ch_b}']
        merged_ch['Adj_CD'] = (
            merged_ch[f'Adjusted_CD_ch{ch_a}'] + merged_ch[f'Adjusted_CD_ch{ch_b}']
        ) / 2.0

        for nerve in sorted(merged_ch['rec_nerve'].unique()):
            sub = merged_ch[merged_ch['rec_nerve'] == nerve]
            first_for_nerve = nerve.lower() not in nerve_legend_shown
            if first_for_nerve:
                nerve_legend_shown.add(nerve.lower())

            hover = [
                f"Exp: {r['expdate']}<br>"
                f"Nerve: {nerve_label(nerve)}<br>"
                f"Pol: {int(r['dataset_stim_polarity'])}<br>"
                f"Avg Adj CD: {r['Adj_CD']:.2f} mm<br>"
                f"PW: {r['stim_PW']:.3f} ms | Amp: {r['max_stim_amplitude_mA']:.2f} mA<br>"
                f"Ch{ch_a}={r[f'{metric_col}_ch{ch_a}']:.4f}, Ch{ch_b}={r[f'{metric_col}_ch{ch_b}']:.4f}<br>"
                f"Diff: {r['Diff']:.4f} {y_unit}"
                for _, r in sub.iterrows()
            ]

            fig.add_trace(go.Scatter(
                x=sub['Adj_CD'], y=sub['Diff'],
                mode='markers',
                name=nerve_label(nerve),
                legendgroup=nerve.lower(),
                showlegend=first_for_nerve,
                marker=dict(
                    color=nerve_color(nerve), symbol='circle',
                    size=9, opacity=0.8,
                    line=dict(color='black', width=1),
                ),
                hovertext=hover, hoverinfo='text',
            ), row=1, col=col_idx)

    # ── y=0 reference lines ────────────────────────────────────────────────
    for col in range(1, 5):
        fig.add_hline(y=0, line_dash='dash', line_color='gray', line_width=1.5, row=1, col=col)

    # ── Axis labels ────────────────────────────────────────────────────────
    for c in range(1, 5):
        fig.update_xaxes(title_text='Adj. CD (mm)', row=1, col=c)
    fig.update_yaxes(title_text=f'Difference ({y_unit})', row=1, col=1)

    fig.update_layout(
        title=f'Paired Differences — {metric_label} vs Adjusted Conduction Distance',
        height=540, width=1700,
        template='plotly_white',
        hovermode='closest',
        font=dict(size=12),
        showlegend=True,
        legend=dict(groupclick='toggleitem'),
    )
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# GENERATE ALL 4 FIGURES
# ──────────────────────────────────────────────────────────────────────────────
print("Generating figures...\n")

fig1 = make_unpaired_figure('peak_amplitude', 'Peak-to-Peak Amplitude (µV)')
fig1.show()
print("Fig 1: Unpaired P2P vs CD")

fig2 = make_unpaired_figure('peak_time', 'Peak Time (ms)')
fig2.show()
print("Fig 2: Unpaired Peak Time vs CD")

fig3 = make_paired_diff_figure('peak_time', 'Peak Time', 'ms')
fig3.show()
print("Fig 3: Paired differences — Peak Time (4 panels)")

fig4 = make_paired_diff_figure('peak_amplitude', 'Peak-to-Peak', 'µV')
fig4.show()
print("Fig 4: Paired differences — P2P (4 panels)")

print("\nDone.")
