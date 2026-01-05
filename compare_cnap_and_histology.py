import numpy as np
import pandas as pd
import os
import yaml
import matplotlib.pyplot as plt

# Helper function to store hist data in a DF - currently stores AUC
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


# Load paths/data
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


# Normalize nerve labels
histology_nerve_map = {
    'sciatic': 'SN',
    'tibial': 'tibial',
    'vagus': 'cVN'
}

hist_auc_df['NerveCanonical'] = hist_auc_df['NerveType'].map(histology_nerve_map)
cnap_meta_df['NerveCanonical'] = cnap_meta_df['rec_nerve']

# Drop anything unmapped (safety)
hist_auc_df = hist_auc_df.dropna(subset=['NerveCanonical'])
cnap_meta_df = cnap_meta_df.dropna(subset=['NerveCanonical'])

# Aggregate CNAP AUC 
avg_auc_df = (
    cnap_meta_df
    .groupby(['expdate', 'NerveCanonical'], as_index=False)
    .agg(
        avg_auc=('auc', 'mean'),
        std_auc=('auc', 'std'),
        n_records=('auc', 'count')
    )
)


# hist -> CNAP mapping
mapping = {
    '39-1': '20220531', #_5-_10
    '39-2': '20220531', #_1-_4
    '40-1': '20220622',
    '41-1': '20220720', #_4-_7
    '41-2': '20220720', #_1-_3
    '43-1': '20221103',
    '44-1': '20221108', #_1-_2
    '44-2': '20221108', #_3-_4
    '45-1': '20221122',
    '46-1': '20221201', #_7-_10
    '46-2': '20221202', #_1-_6
    '47-1': '20230131',
    '48-1': '20230207', #_4-_7
    '48-2': '20230207' #_1-_3

}

hist_auc_df['expdate'] = hist_auc_df['Name'].map(mapping)
hist_auc_df = hist_auc_df.dropna(subset=['expdate'])


# Merge with nerve matching - ensure we only include the right sections from each data
merged_df = pd.merge(
    hist_auc_df,
    avg_auc_df,
    left_on=['expdate', 'NerveCanonical'],
    right_on=['expdate', 'NerveCanonical'],
    how='inner'
)


# sanity checks
print("Merged rows:")
print(merged_df[['Name', 'expdate', 'NerveCanonical']])

print("\nCounts per nerve:")
print(merged_df['NerveCanonical'].value_counts())


# Plot
plt.figure(figsize=(7, 5))
plt.scatter(
    merged_df['hist_auc_fibers_um'],
    merged_df['avg_auc'],
    s=80
)

for _, row in merged_df.iterrows():
    plt.text(
        row['hist_auc_fibers_um'] + 1,
        row['avg_auc'] + 5,
        f"{row['Name']} ({row['NerveCanonical']})",
        fontsize=9
    )

plt.xlabel("Fiber diameter histogram AUC (fibers·μm)")
plt.ylabel("Mean CNAP AUC")
plt.title("Fiber Diameter AUC vs CNAP AUC\n(nerve-matched, temporary mapping)")
plt.grid(True)
plt.tight_layout()
plt.show()