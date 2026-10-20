import pandas as pd
import os
import yaml

# ============================================================
# Load configuration
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = local_config["local_output_path"]

# ============================================================
# Load histology data (axon-level)
# ============================================================
HIST_CSV = os.path.join(LOCAL_OUTPUT_PATH, "all_samples_fiber_diameters.csv")
hist_df = pd.read_csv(HIST_CSV)

# ============================================================
# Compute axon count per histology sample
# ============================================================
axon_counts = (
    hist_df.groupby("Name")
           .agg(axon_count=("SAE_diameter", "count"))
           .reset_index()
)

# ============================================================
# Print results
# ============================================================
print("\n=== Axon Count per Histology Sample ===\n")
for _, row in axon_counts.iterrows():
    print(f"Sample {row['Name']}: {row['axon_count']} axons")
