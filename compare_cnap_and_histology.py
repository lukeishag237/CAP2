import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import yaml
import os

# ============================================================
# Load configuration and CNAP data
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = local_config["local_output_path"]
STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")

cnap_df = pd.read_pickle(STORED_DF_PATH)

# ============================================================
# Load histology data
# ============================================================
RAW_DATA_PATH = local_config.get("raw_data_path")
HIST_CSV = os.path.join(RAW_DATA_PATH, '..', '60x', 'derivatives', 'summary', 'all_samples_fiber_diameters.csv')
hist_df = pd.read_csv(HIST_CSV)
print("Histology columns:", hist_df.columns)

# Map histology nerve names to CNAP nerves
nerve_map = {
    'sciatic': 'SN',
    'tibial': 'tibial',
    'cVN': 'cVN'
}
hist_df['NerveCanonical'] = hist_df['NerveType'].map(nerve_map)

# Map histology sample names to CNAP experiment dates
mapping = {
    '39-1': '20220531', '39-2': '20220531',
    '40-1': '20220622',
    '41-1': '20220720', '41-2': '20220720',
    '43-1': '20221103',
    '44-1': '20221108', '44-2': '20221108',
    '45-1': '20221122',
    '46-1': '20221201', '46-2': '20221202',
    '47-1': '20230131',
    '48-1': '20230207', '48-2': '20230207'
}
hist_df['expdate'] = hist_df['Name'].map(mapping)
hist_df = hist_df.dropna(subset=['expdate', 'NerveCanonical'])

# ============================================================
# Compute histology AUC per nerve/sample (sum of fiber diameters)
# ============================================================
hist_summary = (
    hist_df.groupby(['expdate', 'NerveCanonical'], as_index=False)
    .agg(hist_auc_fibers_um=('SAE_diameter', 'sum'))
)

# ============================================================
# Power-law plots BEFORE filtering for histology
# ============================================================
MAX_DISTANCE_MM = 26
NERVES = ["cVN", "SN", "tibial"]
FIBER_TYPES = ["myel", "unmyel"]
METRICS = {
    "peak_amplitude": "Peak amplitude (mV)",
    "peak_time": "Peak time (ms)",
    "auc": "CNAP AUC"
}
NERVE_COLORS = {
    "cVN": "tab:blue",
    "SN": "tab:green",
    "tibial": "tab:red"
}

# Filter CNAP data for plotting (ch_num_neural = 1, positive distance)
plot_df = cnap_df[
    (cnap_df['rec_nerve'].isin(NERVES)) &
    (cnap_df['ch_num_neural'] == 1) &
    (cnap_df['stim_rec_dist'] > 0)
].copy()

def power_law_fit(x, y):
    x_log = np.log(x)
    y_log = np.log(y)
    X = sm.add_constant(x_log)
    model = sm.OLS(y_log, X).fit()
    k = np.exp(model.params["const"])
    alpha = model.params[x_log.name]
    return k, alpha, model

results = {}
for fiber in FIBER_TYPES:
    results[fiber] = {}
    for nerve in NERVES:
        results[fiber][nerve] = {}
        df_nf = plot_df[
            (plot_df['rec_nerve'] == nerve) &
            (plot_df['target_fiber_types'] == fiber)
        ]
        for metric in METRICS:
            df_m = df_nf[df_nf[metric] > 0]
            if len(df_m) < 5:
                continue
            k, alpha, model = power_law_fit(df_m["stim_rec_dist"], df_m[metric])
            results[fiber][nerve][metric] = {"k": k, "alpha": alpha, "model": model, "df": df_m}
            print(f"{nerve} | {fiber} | {metric}: {metric} = {k:.2f} * distance^{alpha:.2f} (R²={model.rsquared:.3f})")

# Plot power-law fits
d_fit = np.linspace(5, MAX_DISTANCE_MM, 300)
for fiber in FIBER_TYPES:
    for metric, ylabel in METRICS.items():
        plt.figure(figsize=(7,6))
        for nerve in NERVES:
            if metric not in results[fiber][nerve]:
                continue
            df_plot_metric = results[fiber][nerve][metric]['df']
            k = results[fiber][nerve][metric]['k']
            alpha = results[fiber][nerve][metric]['alpha']
            y_fit = k * d_fit ** alpha
            plt.scatter(df_plot_metric["stim_rec_dist"], df_plot_metric[metric],
                        alpha=0.6, color=NERVE_COLORS[nerve], label=nerve)
            plt.plot(d_fit, y_fit, "--", color=NERVE_COLORS[nerve], linewidth=2)
        plt.xlabel("Conduction distance (mm)")
        plt.ylabel(ylabel)
        plt.title(f"{fiber.upper()} fibers: {ylabel} vs distance")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

# ============================================================
# Merge CNAP with histology for regression
# ============================================================
cnap_df_filtered = cnap_df[
    (cnap_df['rec_nerve'].isin(NERVES)) &
    (cnap_df['ch_num_neural'] == 1)
].copy()

cnap_hist_df = pd.merge(
    cnap_df_filtered,
    hist_summary,
    left_on=['expdate', 'rec_nerve'],
    right_on=['expdate', 'NerveCanonical'],
    how='inner'
)

# ============================================================
# Prepare regression dataframe
# ============================================================
df_reg = cnap_hist_df.copy()
df_reg = df_reg[(df_reg['stim_rec_dist'] > 0) & (df_reg['auc'] > 0)].copy()
df_reg["log_dist"] = np.log(df_reg["stim_rec_dist"])
df_reg["log_auc"] = np.log(df_reg["auc"])
finite_mask = np.isfinite(df_reg["log_dist"]) & np.isfinite(df_reg["log_auc"]) & np.isfinite(df_reg['hist_auc_fibers_um'])
df_reg = df_reg.loc[finite_mask].copy()

# Categorical variables
df_reg["fiber"] = df_reg["target_fiber_types"].astype("category")
df_reg["nerve"] = df_reg["rec_nerve"].astype("category")
df_reg["subject"] = df_reg["expdate"].astype("category")

# ============================================================
# LINEAR MODEL (including histology)
# ============================================================
X_lin = df_reg[['log_dist', 'fiber', 'nerve', 'hist_auc_fibers_um']]
X_lin = pd.get_dummies(X_lin, drop_first=True)
X_lin = sm.add_constant(X_lin)
y_lin = df_reg['log_auc']
linear_model = sm.OLS(y_lin, X_lin.astype(float)).fit()
print("\n================ LINEAR MODEL ================\n")
print(linear_model.summary())

# ============================================================
# FIXED EFFECTS MODEL (subject dummies)
# ============================================================
X_fe = df_reg[['log_dist', 'fiber', 'nerve', 'subject', 'hist_auc_fibers_um']]
X_fe = pd.get_dummies(X_fe, drop_first=True)
X_fe = sm.add_constant(X_fe)
y_fe = df_reg['log_auc']
fe_model = sm.OLS(y_fe, X_fe.astype(float)).fit()
print("\n================ FIXED EFFECTS MODEL ================\n")
print(fe_model.summary())

# ============================================================
# MIXED EFFECTS MODEL (random subject intercept)
# ============================================================
me_model = smf.mixedlm(
    "log_auc ~ log_dist + fiber + nerve + hist_auc_fibers_um",
    df_reg,
    groups=df_reg["subject"]
)
me_result = me_model.fit(method="lbfgs")
print("\n================ MIXED EFFECTS MODEL ================\n")
print(me_result.summary())
