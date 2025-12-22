import numpy as np
import pandas as pd
import os
import yaml
import matplotlib.pyplot as plt

def compute_histogram_auc_per_nerve(csv_file: str, bin_width: float = 0.07) -> pd.DataFrame:
    """
    Compute histogram AUC for each nerve using the same binning
    logic as plot_diameter_histograms_grid().
    
    Returns:
        DataFrame with one row per nerve.
    """
    df = pd.read_csv(csv_file)

    results = []

    for (nerve_type, name), group in df.groupby(['NerveType', 'Name']):
        diameters = group['MinFeret'].values

        if len(diameters) == 0:
            continue

        max_diameter = diameters.max()
        bin_edges = np.arange(0, max_diameter + bin_width, bin_width)

        counts, _ = np.histogram(diameters, bins=bin_edges)

        # Histogram AUC (counts × bin width)
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


CONFIG_FILE = "local_config.yml"

try:
    with open(CONFIG_FILE, 'r') as f:
        # Use yaml.safe_load() for security
        local_config = yaml.safe_load(f) 
except FileNotFoundError:
    raise FileNotFoundError(f"Config file '{CONFIG_FILE}' not found. "
                          f"Did you copy 'config.yml.example' to 'config.yml'?")

RAW_DATA_PATH = local_config.get("raw_data_path")

csv_file = os.path.join(RAW_DATA_PATH, '..', '60x', 'derivatives', 'summary', 'all_samples_fiber_diameters.csv')

# Only takes 6 seconds so decided not to use a df intermediate
hist_auc_df = compute_histogram_auc_per_nerve(csv_file= csv_file)

LOCAL_OUTPUT_PATH = local_config.get("local_output_path")
STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH,"stored_df.pkl")

cnap_meta_df = pd.read_pickle(STORED_DF_PATH)
avg_auc_df = (cnap_meta_df.groupby('expdate', as_index=False).agg(avg_auc=('auc', 'mean'),std_auc=('auc', 'std'),n_records=('auc', 'count')).sort_values('expdate'))

# print(hist_auc_df)
# print(avg_auc_df)

# Assume both df's order is correct: plot avg AUC of fiber diams vs. cnaps auc
mapping = {
    '39-2': '20220511',
    '41-2': '20220518',
    '44-1': '20220525',
    '46-2': '20220531',
    '48-2': '20220616',
    '39-1': '20220622',
    '41-1': '20220720',
    '44-2': '20221004',
    '46-1': '20221103',
    '48-1': '20221108',
    '40-1': '20221122',
    '43-1': '20221201',
    '45-1': '20230131',
    '47-1': '20230207'
}

# Map expdates to histology df
hist_auc_df['expdate'] = hist_auc_df['Name'].map(mapping)

# Merge with CNAP avg AUC
merged_df = pd.merge(
    hist_auc_df,
    avg_auc_df[['expdate', 'avg_auc']],
    on='expdate',
    how='inner'
)

# --- Plot ---
plt.figure(figsize=(7,5))
plt.scatter(merged_df['hist_auc_fibers_um'], merged_df['avg_auc'], s=80)
for i, row in merged_df.iterrows():
    plt.text(row['hist_auc_fibers_um']+1, row['avg_auc']+5, row['Name'], fontsize=9)
plt.xlabel("Fiber diameter histogram AUC (fibers·μm)")
plt.ylabel("Mean CNAP AUC")
plt.title("Fiber Diameter AUC vs CNAP AUC (temporary mapping)")
plt.grid(True)
plt.tight_layout()
plt.show()