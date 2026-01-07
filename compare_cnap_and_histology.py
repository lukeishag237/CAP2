import numpy as np
import pandas as pd
import os
import yaml
import matplotlib.pyplot as plt


# ============================================================
# Helper: compute histology histogram AUC per nerve / sample
# ============================================================
def compute_histogram_auc_per_nerve(csv_file: str, bin_width: float = 0.07) -> pd.DataFrame:
    df = pd.read_csv(csv_file)
    results = []

    for (nerve_type, name), group in df.groupby(['NerveType', 'Name']):
        diameters = group['MinFeret'].values
        if len(diameters) == 0:
            continue

        max_diameter = diameters.max()
        bin_edges = np.arange(0, max_diameter + bin_width, bin_width)
        counts, _ = np.histogram(diameters, bins=bin_edges)

        auc = np.sum(counts) * bin_width

        results.append({
            'NerveType': nerve_type,
            'Name': name,
            'n_fibers': len(diameters),
            'bin_width_um': bin_width,
            'hist_auc_fibers_um': auc,
            'fiber_count_check': np.sum(counts)
        })

    return pd.DataFrame(results)


# ============================================================
# Load paths and data
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, 'r') as f:
    local_config = yaml.safe_load(f)

RAW_DATA_PATH = local_config["raw_data_path"]
LOCAL_OUTPUT_PATH = local_config["local_output_path"]

csv_file = os.path.join(
    RAW_DATA_PATH, '..', '60x', 'derivatives', 'summary',
    'all_samples_fiber_diameters.csv'
)

STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")

hist_auc_df = compute_histogram_auc_per_nerve(csv_file)
cnap_meta_df = pd.read_pickle(STORED_DF_PATH)


# ============================================================
# Normalize nerve labels
# ============================================================
histology_nerve_map = {
    'sciatic': 'SN',
    'tibial': 'tibial',
    'vagus': 'cVN'
}

hist_auc_df['NerveCanonical'] = hist_auc_df['NerveType'].map(histology_nerve_map)
cnap_meta_df['NerveCanonical'] = cnap_meta_df['rec_nerve']

hist_auc_df = hist_auc_df.dropna(subset=['NerveCanonical'])
cnap_meta_df = cnap_meta_df.dropna(subset=['NerveCanonical'])


# ============================================================
# Aggregate CNAP AUC (for averaged plot)
# ============================================================
avg_auc_df = (
    cnap_meta_df
    .groupby(['expdate', 'NerveCanonical'], as_index=False)
    .agg(
        avg_auc=('auc', 'mean'),
        std_auc=('auc', 'std'),
        n_records=('auc', 'count')
    )
)


# ============================================================
# Histology → CNAP date mapping
# ============================================================
mapping = {
    '39-1': '20220531',
    '39-2': '20220531',
    '40-1': '20220622',
    '41-1': '20220720',
    '41-2': '20220720',
    '43-1': '20221103',
    '44-1': '20221108',
    '44-2': '20221108',
    '45-1': '20221122',
    '46-1': '20221201',
    '46-2': '20221202',
    '47-1': '20230131',
    '48-1': '20230207',
    '48-2': '20230207'
}

hist_auc_df['expdate'] = hist_auc_df['Name'].map(mapping)
hist_auc_df = hist_auc_df.dropna(subset=['expdate'])


# ============================================================
# Merge for averaged plot
# ============================================================
merged_avg_df = pd.merge(
    hist_auc_df,
    avg_auc_df,
    on=['expdate', 'NerveCanonical'],
    how='inner'
)


# ============================================================
# Plot 1: histology AUC vs mean CNAP AUC
# ============================================================
plt.figure(figsize=(7, 5))

plt.scatter(
    merged_avg_df['hist_auc_fibers_um'],
    merged_avg_df['avg_auc'],
    s=80
)

for _, row in merged_avg_df.iterrows():
    plt.text(
        row['hist_auc_fibers_um'] + 1,
        row['avg_auc'] + 5,
        f"{row['Name']} ({row['NerveCanonical']})",
        fontsize=9
    )

plt.xlabel("Fiber diameter histogram AUC (fibers·μm)")
plt.ylabel("Mean CNAP AUC")
plt.title("Fiber Diameter AUC vs CNAP AUC\n(histology-averaged, nerve-matched)")
plt.grid(True)
plt.tight_layout()
plt.show()


# ============================================================
# Plot 2: per-record CNAP vs histology AUC (ch_num_neural == 1)
# ============================================================
valid_hist_pairs = hist_auc_df[['expdate', 'NerveCanonical']].drop_duplicates()

# Filter CNAP rows:
#   1) must have histology
#   2) must be neural channel 1
cnap_filtered = pd.merge(
    cnap_meta_df,
    valid_hist_pairs,
    on=['expdate', 'NerveCanonical'],
    how='inner'
)

cnap_filtered = cnap_filtered[cnap_filtered['ch_num_neural'] == 1]

# Attach histology AUC
cnap_with_hist = pd.merge(
    cnap_filtered,
    hist_auc_df[['expdate', 'NerveCanonical', 'hist_auc_fibers_um']],
    on=['expdate', 'NerveCanonical'],
    how='inner'
)

plt.figure(figsize=(7, 5))

plt.scatter(
    cnap_with_hist['hist_auc_fibers_um'],
    cnap_with_hist['auc'],
    s=60,
    alpha=0.7
)

plt.xlabel("Fiber diameter histogram AUC (fibers·μm)")
plt.ylabel("CNAP AUC (per recording)")
plt.title("Fiber Diameter AUC vs CNAP AUC\n(per-record, ch_num_neural = 1)")
plt.grid(True)
plt.tight_layout()
plt.show()
