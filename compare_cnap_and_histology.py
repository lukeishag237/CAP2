import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import yaml
import os

# ============================================================
# Load data
# ============================================================
CONFIG_FILE = "local_config.yml"

with open(CONFIG_FILE, "r") as f:
    local_config = yaml.safe_load(f)

LOCAL_OUTPUT_PATH = local_config["local_output_path"]
STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH, "stored_df.pkl")

df = pd.read_pickle(STORED_DF_PATH)

# ============================================================
# Analysis parameters
# ============================================================
MAX_DISTANCE_MM = 26  # match Fig 8 short-distance regime

NERES = ["cVN", "SN", "tibial"]
FIBER_TYPES = ["myel", "unmyel"]

METRICS = {
    "peak_amplitude": "Peak amplitude (mV)",
    "peak_time": "Peak time (ms)",
    "auc": "CNAP AUC"
}

NERve_COLORS = {
    "cVN": "tab:blue",
    "SN": "tab:green",
    "tibial": "tab:red"
}

FIBER_MARKERS = {
    "myel": "o",
    "unmyel": "s"
}

# ============================================================
# Filter base dataframe
# ============================================================
df = df[
    (df["ch_num_neural"] == 1) &
    (df["stim_rec_dist"] > 0) &
    (df["stim_rec_dist"] < MAX_DISTANCE_MM) &
    (df["rec_nerve"].isin(NERES)) &
    (df["target_fiber_types"].isin(FIBER_TYPES))
].copy()

# ============================================================
# Power-law regression helper
# ============================================================
def power_law_fit(x, y):
    """
    Fits y = k * x^alpha via log-log regression
    Returns k, alpha, fitted model
    """
    x_log = np.log(x)
    y_log = np.log(y)

    X = sm.add_constant(x_log)
    model = sm.OLS(y_log, X).fit()

    k = np.exp(model.params["const"])
    alpha = model.params[x_log.name]

    return k, alpha, model

# ============================================================
# Fit models
# ============================================================
results = {}

for nerve in NERES:
    results[nerve] = {}

    for fiber in FIBER_TYPES:
        results[nerve][fiber] = {}

        df_nf = df[
            (df["rec_nerve"] == nerve) &
            (df["target_fiber_types"] == fiber)
        ]

        for metric in METRICS:
            df_m = df_nf[df_nf[metric] > 0]

            if len(df_m) < 5:
                continue

            k, alpha, model = power_law_fit(
                df_m["stim_rec_dist"],
                df_m[metric]
            )

            results[nerve][fiber][metric] = {
                "k": k,
                "alpha": alpha,
                "model": model,
                "df": df_m
            }

            print(
                f"{nerve} | {fiber} | {metric}: "
                f"{metric} = {k:.2f} * distance^{alpha:.2f} "
                f"(R²={model.rsquared:.3f})"
            )

# ============================================================
# Plotting
# ============================================================
d_fit = np.linspace(5, MAX_DISTANCE_MM, 300)

for metric, ylabel in METRICS.items():
    plt.figure(figsize=(7, 6))

    for nerve in NERES:
        for fiber in FIBER_TYPES:
            if metric not in results[nerve][fiber]:
                continue

            df_plot = results[nerve][fiber][metric]["df"]
            k = results[nerve][fiber][metric]["k"]
            alpha = results[nerve][fiber][metric]["alpha"]

            y_fit = k * d_fit ** alpha

            plt.scatter(
                df_plot["stim_rec_dist"],
                df_plot[metric],
                alpha=0.6,
                color=NERve_COLORS[nerve],
                marker=FIBER_MARKERS[fiber],
                label=f"{nerve}-{fiber}"
            )

            plt.plot(
                d_fit,
                y_fit,
                color=NERve_COLORS[nerve],
                linestyle="--",
                linewidth=2
            )

    plt.xlabel("Conduction distance (mm)")
    plt.ylabel(ylabel)
    plt.title(f"Effect of conduction distance on {ylabel}")
    plt.grid(True)

    # Deduplicate legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), fontsize=8)

    plt.tight_layout()
    plt.show()
