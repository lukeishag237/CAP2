CAP2 — CAP / CNAP preprocessing and analysis

Overview
----
This repository contains MATLAB scripts to preprocess compound action potential (CAP) recordings and a Python analysis script to load, extract features from, and visualize the preprocessed results.

Key features
----
- MATLAB preprocessing pipeline and example analysis scripts for collecting and organizing CAP data.
- A Python analyzer (`cnap_analyzer.py`) that loads preprocessed MATLAB `.mat` files, computes features (peak amplitude, peak latency, AUC), validates signals visually, and produces summary plots/metrics.
- Configuration via `local_config.yml` to point the code at local raw and preprocessed data locations.

Repository structure (important files)
----
- `AnalysisScript_CAPs_Batch.m`, `AnalysisScript_CAPs_Wrapper.m`, `AnalysisScript_CAPs.m` — MATLAB scripts used to preprocess and batch-process CAP datasets and organize outputs.
- `cnap_analyzer.py` — Python analysis and plotting tool; loads preprocessed `.mat` files, computes CAP/CNAP features, and creates figures and summary metrics.
- `local_config.example.yml` — Example YAML config; copy to `local_config.yml` and edit paths for your environment.
- `biomed-base.yml` — Conda environment YAML listing recommended dependencies for running the Python analysis.
- `external_dependencies/adinstruments_sdk_matlab-master/` — Third-party MATLAB SDK for reading LabChart/PowerLab files (used by preprocessing scripts).
- `excluded_datasets.csv` — CSV of datasets/trials explicitly excluded from analysis.
- `local_output/` — Default folder for generated figures and metrics

Quick start
----
1. Create a local config
	- Copy `local_config.example.yml` to `local_config.yml` and set absolute paths for `raw_data_path`, `preprocessed_CAP_path`, and (optionally) `local_output_path`.

2. Install dependencies (recommended using conda/mamba)
	- Create the environment from `biomed-base.yml`:

```powershell
mamba env create --name biomed-base --file biomed-base.yml
conda activate biomed-base
```

3. Preprocess (MATLAB)
	- Run the MATLAB preprocessing scripts to generate the preprocessed `.mat` files expected by the Python analyzer. The MATLAB scripts use the `external_dependencies/adinstruments_sdk_matlab-master/` package to read LabChart files.
    - NOTE: This part only works in Windows because the LabChart SDK only runs in Windows

4. Run Python analysis
	- Update `local_config.yml` as above, then run the `cnap_analyzer.py` script from the repository root:

```powershell
python cnap_analyzer.py
```

This will load the preprocessed `.mat` files (from `preprocessed_CAP_path`), compute features, populate `processor.data_df`, and generate figures and metrics in `local_output/`.

Configuration details
----
- `local_config.yml`:
  - `raw_data_path`: path to unprocessed raw recordings and metadata spreadsheet.
  - `preprocessed_CAP_path`: path where MATLAB preprocessing writes `.mat` files.
  - `local_output_path`: optional path for saving plots and metrics (defaults to `local_output`).

Notable code points
----
- `cnap_analyzer.CNAPAnalyzer` — Encapsulates signal-windowing logic, peak detection, and AUC calculation.
- `cnap_analyzer.CNAPDataProcessor` — Loads metadata, iterates preprocessed files, computes features, and generates summary plots/grids.
- MATLAB scripts under the repo root handle the heavy lifting of reading raw LabChart data and producing the `.mat` files consumed by the Python analyzer.

External dependencies
----
- The MATLAB ADInstruments SDK (included under `external_dependencies/`) is required by the MATLAB preprocessing scripts.
- Python dependencies are listed in `biomed-base.yml` (numpy, scipy, pandas, matplotlib, seaborn, etc.).

Troubleshooting and notes
----
- If `cnap_analyzer.py` raises a `FileNotFoundError` for `local_config.yml`, ensure you copied the example config and set absolute paths.
- Paths in `local_config.yml` should point to existing directories — both the raw data and preprocessed-cap folders are asserted at script start.
- The repository assumes MATLAB preprocessing has produced `.mat` files with expected variable names (e.g., `t_window`, `ch_neural_avg`). If you modify preprocessing outputs, update `cnap_analyzer.py` accordingly.
