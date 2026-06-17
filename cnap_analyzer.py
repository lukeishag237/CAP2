import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy import integrate
from pathlib import Path
import glob
from typing import Dict, List, Tuple, Optional
import yaml

class CNAPAnalyzer:
    """
    A class with methods to calculate CAP features (e.g., positive peak, area
    under the curve (AUC)) based on specified CV bounds and artifact bounds, and
    a method to plot each CAP with bounds and positive peak overlaid as gut check
    """
    
    def __init__(self, 
                 sampling_rate: float = np.nan, # sampling rate in Hz; set to nan until set per recording
                 artifact_window_ms: float = 0.16, # 2.0,
                 min_cv_ms: float = 5, #0.2,  # minimum conduction velocity in m/s; used to set max bound
                 distance_mm: float = np.nan,  # distance between stimulus and recording in mm; set to nan until set per recording
                 metadata_df = None # loaded data from the xlsx file - we will use this to calculate len of stim artifact
                ):
        """
        Initialize the CNAP analyzer.
        
        Args:
            sampling_rate: Sampling rate in Hz
            artifact_window_ms: Time window after stimulus to exclude (artifact duration)
            min_cv_ms: Minimum conduction velocity to consider (m/s)
            max_cv_ms: Maximum conduction velocity to consider (m/s)
            distance_mm: Distance between stimulus and recording site (mm)
        """
        # Set sampling rate to not a number initially since it should not be used until set per recording
        self.sampling_rate = sampling_rate
        self.artifact_window_ms = artifact_window_ms
        self.min_cv = min_cv_ms
        self.distance = distance_mm
        self.metadata_df = metadata_df # we will set this in the dataprocessor.plot_data method
        self.idx = np.nan # keep track of which row in metadata df we are in
        

    #NEW: Combined with peak finding to make sure we don't go too far
    def find_signal_bounds_and_peak(self, signal: np.ndarray, times: np.ndarray) -> Tuple[int, int, Dict[str, float], ]:
        """
        Find the appropriate start and end indices for signal analysis.
        
        Args:
            signal: CNAP signal array
            times: Time points array (ms)
            
        Returns:
            Tuple of (start_idx, end_idx)
        """
        # Find artifact end (start of analysis window)
        polarity = 1
        if self.metadata_df['phase_type'][self.idx] == "biphasic":
            polarity = 2
        
        # Find analysis window end based on conduction velocity. Added an if statement to implement unmyel fibers
        if self.metadata_df["target_fiber_types"][self.idx] == 'myel':
            analysis_end_idx = np.where(times >= (self.distance / self.min_cv))[0][0]
            artifact_end_ms = float(self.metadata_df['stim_PW'][self.idx]) * polarity + 0.17 # increased to .2 from .16 due to peak time error
        else: #unmyel
            analysis_end_idx = np.where(times >= (self.distance / .5))[0][0] # what should this value be?
            artifact_end_ms = float(self.metadata_df['stim_PW'][self.idx]) * polarity + 2 # is this wrong?
        
        # Calculate important indexes to return
        artifact_end_idx = np.where(times >= artifact_end_ms)[0][0]
        peak_idx = artifact_end_idx + np.argmax(signal[artifact_end_idx:analysis_end_idx])

        # New! if CAP is DECREASING at artifact end (analysis start) - we must've went too far
        # Backtrack to peak and see if it is higher than our current peak  - reset artifact idx if it is
        if 1 < artifact_end_idx < len(signal) - 1:

            slope = signal[artifact_end_idx] - signal[artifact_end_idx - 1]
            if slope < -0.1:   # much smaller threshold needed
                i = artifact_end_idx - 1
                # Walk backward until slope flips or start is reached
                while i > 1:
                    prev_slope = signal[i] - signal[i - 1]
                    next_slope = signal[i + 1] - signal[i]

                    # Local maximum condition
                    if prev_slope > 0 and next_slope <= 0:
                        # Only update if this peak is higher than current peak
                        if signal[i] > signal[peak_idx]:
                            artifact_end_idx = i-10
                            peak_idx = i
                        break

                    i -= 1

        return artifact_end_idx, analysis_end_idx, {'peak_time': times[peak_idx],'peak_amplitude': signal[peak_idx]}
    
    # Rolled this functionality into the find_signal_bounds function - unused now
    # def find_peak(self, signal: np.ndarray, times: np.ndarray) -> Dict[str, float]:
    #     """
    #     Find the maximum peak in the neural component of the signal.     
    #     Args:
    #         signal: CNAP signal array
    #         times: Time points array (ms) 
    #     Returns:
    #         Dictionary containing peak time and amplitude
    #     """
    #     start_idx, end_idx = self.find_signal_bounds(signal, times)
    #     peak_idx = start_idx + np.argmax(signal[start_idx:end_idx])
    #     return {'peak_time': times[peak_idx],'peak_amplitude': signal[peak_idx]}
    
    def calculate_auc(self, 
                     signal: np.ndarray, 
                     times: np.ndarray
                    ) -> float:
        """
        Calculate area under the curve for the neural component.
        
        Args:
            signal: CNAP signal array
            times: Time points array (ms)
            method: Integration method ('trapz' or 'simpson')
            
        Returns:
            Area under the curve
        """
        start_idx, end_idx, peak_info = self.find_signal_bounds_and_peak(signal, times)
        
        # auc = integrate.trapz(np.abs(signal[start_idx:end_idx]),
        #                         times[start_idx:end_idx])
        auc = integrate.trapezoid(np.abs(signal[start_idx:end_idx]),
                                    times[start_idx:end_idx])
            
        return auc
    
    def plot_gut_check_bounds_and_peak(self, 
                     signal: np.ndarray, 
                     times: np.ndarray,
                     show_bounds: bool = True
                    ) -> None:
        """
        Plot the signal with analysis annotations.
        
        Args:
            signal: CNAP signal array
            times: Time points array (ms)
            show_bounds: Whether to show analysis bounds
        """
        start_idx, end_idx, peak_info = self.find_signal_bounds_and_peak(signal, times)
        
        plt.figure(figsize=(10, 6))
        plt.plot(times, signal, 'b-', label='CNAP')
        
        if show_bounds:
            plt.axvline(x=times[start_idx], color='r', linestyle='--', 
                       label='Analysis Start')
            plt.axvline(x=times[end_idx], color='g', linestyle='--', 
                       label='Analysis End')
        
        plt.plot(peak_info['peak_time'], peak_info['peak_amplitude'], 'ro', 
                label='Peak')
        

        buffer_time_to_plot = (times[end_idx]-times[start_idx])*0.05
        plot_start_idx = np.where(times == 0)[0][0] - 10
        # plt.xlim(times[start_idx] - buffer_time_to_plot, times[end_idx] + buffer_time_to_plot)
        plt.xlim(times[plot_start_idx], times[end_idx] + buffer_time_to_plot)
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        id = self.metadata_df["file_prefix_outputs_values"][self.idx]
        plt.title(f'CNAP Analysis - {id}')
        plt.legend()
        plt.grid(True)
        plt.show()

class CNAPDataProcessor:
    """
    A class with methods to load the preprocessed CAP signals from .mat
    files, to process CNAP data, and visualize results.
    Uses CNAPAnalyzer to get features.
    """
    def __init__(self, metadata_file: str, local_output_path: str, loaded_df = None):
        """
        Initialize with path to metadata Excel file.
        
        Args:
            metadata_file (str): Path to Excel file containing CNAP metadata
        """
        self.metadata_file = metadata_file
        self.data_df = loaded_df
        # Initialize the CNAP analyzer with default parameters
        self.analyzer = CNAPAnalyzer()
        self.output_path = local_output_path
        
    def load_mat_file(self, mat_file: str, flag_plot: bool = False) -> Dict:
        """
        Load max CAP from specified .mat file & use CNAPAnalyzer to get features.
        Return max CAP along with features.
        
        Args:
            mat_file (str): Path to .mat file
        
        Returns:
            dict: Processed data from .mat file including CNAP analysis
        """
        # Load MATLAB data
        mat_data = loadmat(mat_file)
        
        # Extract time window and signal data
        t_window = mat_data['t_window'].flatten()
        ch_neural_avg = mat_data['ch_neural_avg']
        
        # Update analyzer parameters for this recording
        self.analyzer.distance = float(mat_data['stim_rec_dist'][0][0])
        # Set the sampling rate based on the difference between t_window points
        self.analyzer.sampling_rate = 1 / np.mean(np.diff(t_window)) * 1000  # Convert to Hz from ms
        # Round sampling rate to nearest integer for consistency
        self.analyzer.sampling_rate = int(round(self.analyzer.sampling_rate))

        # Initialize storage for analysis results
        peak_info = {}
        auc_values = {}
        
        # Process each channel if multiple channels exist
        if len(ch_neural_avg.shape) > 1:
            max_amp_idx = np.argmax(mat_data['amp_values_actual'].flatten())
            signal = ch_neural_avg[:, max_amp_idx]
        else:
            signal = ch_neural_avg
        
        # Perform CNAP analysis
        analysis_start,analysis_end, peak_info = self.analyzer.find_signal_bounds_and_peak(signal, t_window)
        auc = self.analyzer.calculate_auc(signal, t_window)

        # Plot individual CAP with positive peak and bounds (optional)
        if flag_plot:
            self.analyzer.plot_gut_check_bounds_and_peak(signal, t_window)

        # Build return dictionary with original and new analysis data
        return {
            'expdate': str(mat_data['expdate'][0]),
            't_window': t_window,
            'ch_neural_avg': signal,
            'ch_num_neural': int(mat_data['ch_num_neural'][0][0]),
            'stim_rec_dist': float(mat_data['stim_rec_dist'][0][0]),
            # Add new CNAP analysis results
            'peak_time': peak_info['peak_time'],
            'peak_amplitude': peak_info['peak_amplitude'],
            'auc': auc
        }

    def process_data(self, preprocessed_CAP_path: str, flag_plot: bool = False) -> None:
        """Load all CAPs and store CAPs and features in a DataFrame."""
        metadata_df = pd.read_excel(self.metadata_file)
        self.analyzer.metadata_df = metadata_df
        data_list = []
        
        for idx, row in metadata_df.iterrows():
            self.analyzer.idx = idx
            file_prefix = row['file_prefix_outputs_values']
            mat_files = glob.glob(str(Path(preprocessed_CAP_path) / file_prefix / '*.mat'))
            
            # print warning if no .mat files found for this row
            if len(mat_files) == 0:
                print(f'Warning: No .mat files found for file prefix {file_prefix}')
                continue

            for mat_file in mat_files:
                data_dict = self.load_mat_file(mat_file,flag_plot)
                data_dict.update({
                    'dataset_stim_polarity': row['dataset_stim_polarity'],
                    'rec_nerve': row['rec_nerve'],
                    'target_fiber_types': row['target_fiber_types']
                })
                data_list.append(data_dict)
        
        self.data_df = pd.DataFrame(data_list)

    # Plot all CAPs with specified experiment date, fiber type, and nerve type
    def plot_data(self, exp_date: str='20220511', fiber_type: str='myel', 
                 nerve: str='SN', show_analysis: bool=True) -> None:
        """
        Create CNAP visualization.
        Similar to MATLAB's gramm plotting but using seaborn's FacetGrid.
        
        Args:
            exp_date (str): Experiment date to plot
            fiber_type (str): Type of fiber to plot
            nerve (str): Nerve type to plot
        """
        # Subset of data to plot must match specified experiment date, target
        # fiber type, and target nerve type
        subset = self.data_df[
            (self.data_df['expdate'].str.contains(exp_date)) & 
            (self.data_df['target_fiber_types'] == fiber_type) & 
            (self.data_df['rec_nerve'] == nerve)
        ]
        
        # Plot the data as follows: Rows are the polarity, columns are the
        # contact number, and line color is conduction distance
        g = sns.FacetGrid(
            data=subset,
            row='dataset_stim_polarity',
            col='ch_num_neural',
            height=4,
            aspect=1.5                                             
        )
        
        def plot_lines(data, **kwargs):
            """Plot lines for each channel with analysis annotations"""
            min_dist = data['stim_rec_dist'].min()
            max_dist = data['stim_rec_dist'].max()
            norm = lambda x: (x - min_dist) / (max_dist - min_dist) if max_dist != min_dist else 0.5
            
            for _, row in data.iterrows():
                plt.plot(
                    row['t_window'],
                    row['ch_neural_avg'],
                    label=f"Dist: {row['stim_rec_dist']}mm",
                    color=plt.cm.viridis(norm(row['stim_rec_dist']))
                )
                
                if show_analysis:
                    # Add peak marker
                    plt.plot(row['peak_time'], row['peak_amplitude'], 'ro', 
                            markersize=5)
                    # Add AUC annotation
                    plt.text(0.7, 0.9, f"AUC: {row['auc']:.2f}", 
                            transform=plt.gca().transAxes)
        

        # Apply plotting function to each subplot
        g.map_dataframe(plot_lines)
        
        # Customize plot
        g.set_axis_labels('time (ms)', 'V$_{pk-pk}$ ($\mu$V)')  
        g.set_titles(col_template="ch #{col_name}", row_template="polarity {row_name}")
        if (fiber_type == 'myel'):
            g.axes[0,0].set_xlim([-0.1, 1.5])
        elif (fiber_type == 'unmyel'):
            g.axes[0,0].set_xlim([-0.1,20])
        else: 
            raise ValueError("Unknown Fiber Type!")
        
        norm = plt.Normalize(subset['stim_rec_dist'].min(), subset['stim_rec_dist'].max())
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        # deprecated; use new form of: cbar_ax = g.fig.add_axes([.92, .3, .02, .4])
        cbar_ax = g.figure.add_axes([.92, .3, .02, .4])
        g.figure.colorbar(sm, cax=cbar_ax, label='stim rec distance (mm)')
        
        # Write the expdate, nerve, and fiber type
        print(f'Exp Date: {exp_date}, Nerve: {nerve}, Fiber Type: {fiber_type}')



        # plt.tight_layout()
        plt.show()
    
    # DEfine a function to plot the peak amplitude vs. condcution distance for the specified exp date, fiber type, and nerve.
    # Allow the peak latency or AUC to be plotted instead of the peak amplitude.
    # def plot_peak_amplitude_vs_distance(self, exp_date: str, fiber_type: str, nerve: str) -> None:
    def plot_peak_amplitude_vs_distance(self, exp_date: str, fiber_type: str, nerve: str, y_var: str = 'peak_amplitude') -> None:
        """Plot peak amplitude vs. conduction distance for the specified exp date, fiber type, and nerve."""
        
        # Filter data for the specified exp date, fiber type, and nerve.
        subset = self.data_df[(
            self.data_df['expdate'].str.contains(exp_date)) & (
            self.data_df['target_fiber_types'] == fiber_type) & (
            self.data_df['rec_nerve'] == nerve)
        ]

        # Plot peak amplitude vs. conduction distance.
        # Use different colors for the dataset_stim_polarity and use different markers  for the contact number.
        # Use colors that are easy to distinguish. The default colors are not easy to distinguish.
        plt.figure(figsize=(10, 6))
        # sns.scatterplot(data=subset, x='stim_rec_dist', y='peak_amplitude', hue='dataset_stim_polarity', style='ch_num_neural', palette='Set1', markers=True, alpha=0.8, s=75)
        sns.scatterplot(data=subset, x='stim_rec_dist', y=y_var, hue='dataset_stim_polarity', style='ch_num_neural', palette='Set1', markers=True, alpha=0.8, s=75)
        # DO the above with a swarmplot
        # sns.swarmplot(data=subset, x='stim_rec_dist', y='peak_amplitude', hue='dataset_stim_polarity', style='ch_num_neural', palette='Set1', markers=True)
        plt.xlabel('Conduction Distance (mm)')

        # plt.ylabel('Peak Amplitude ($\mu$V)')
        plt.ylabel(y_var)
        # plt.title(f'Peak Amplitude vs. Conduction Distance for {nerve} {fiber_type}')
        # plt.legend(title='Stim Polarity')
        plt.grid(True)
        plt.show()

    def plot_metric_grid(self, metric: str, nerve: str, fiber_type: str, ncols: int = 5, save_fig: bool = False) -> None:
        """
        Create a grid of plots for a given metric (peak_amplitude, peak_time, or auc)
        for a specific nerve type, with a fixed number of columns.
        
        Args:
            metric (str): The metric to plot ('peak_amplitude', 'peak_time', or 'auc')
            nerve (str): The nerve type to plot
            ncols (int): Number of columns in the grid
            fiber_type (str): myel or unmyel
            save_fig (bool): Whether to save the figure
        """
        # Get dates for this nerve
        dates = sorted(self.data_df[self.data_df['rec_nerve'] == nerve]['expdate'].unique())
        
        # Calculate number of rows needed
        nrows = (len(dates) + ncols - 1) // ncols  # Ceiling division
        
        # Create figure and axes grid with extra space for title and legend
        # Space in in inches, so we need to convert from inches to pixels (96 dpi) to get the right size.
        fig = plt.figure(figsize=(15, 5),dpi=96) 
        gs = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.1, right=0.82, left=0.1)  # Reduce width to make room for legend
        axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(ncols)] for i in range(nrows)])
        
        # Dictionary for metric labels
        metric_labels = {
            'peak_amplitude': 'Peak Amplitude (µV)',
            'peak_time': 'Peak Time (ms)',
            'auc': 'AUC'
        }
        
        # Find global x-axis limits across all plots (filter by nerve and myelinated fibers)
        nerve_data = self.data_df[
            (self.data_df['rec_nerve'] == nerve) &
            (self.data_df['expdate'].isin(dates)) &
            (self.data_df['target_fiber_types'] == fiber_type)
        ]

        # Commented this out because it is no longer true!
        # Print warning that only myelinated fibers are being plotted for this nerve
        # print(f"Warning: Only myelinated fibers are being plotted for {nerve}.")

        if not nerve_data.empty:
            xmin = nerve_data['stim_rec_dist'].min()
            xmax = nerve_data['stim_rec_dist'].max()
            # Add 10% padding to x limits
            x_padding = (xmax - xmin) * 0.1
            global_xlims = (xmin - x_padding, xmax + x_padding)
        else:
            global_xlims = (0, 1)  # Default limits if no data

        # Find global y-axis limits for each row
        row_ylims = []
        for row in range(nrows):
            start_idx = row * ncols
            end_idx = min((row + 1) * ncols, len(dates))
            row_dates = dates[start_idx:end_idx]
            
            row_data = self.data_df[
                self.data_df['expdate'].isin(row_dates)
            ][metric]
            
            if not row_data.empty:
                ymin = row_data.min()
                ymax = row_data.max()
                # Add 10% padding to the limits
                padding = (ymax - ymin) * 0.1
                row_ylims.append((ymin - padding, ymax + padding))
            else:
                row_ylims.append((0, 1))  # Default limits if no data
        
        # Plot data for each date
        for idx, date in enumerate(dates):
            row = idx // ncols
            col = idx % ncols
            ax = axes[row, col]
            
            # Filter data for current date
            subset = self.data_df[self.data_df['expdate'] == date]
            
            if not subset.empty:
                sns.scatterplot(
                    data=subset,          # DataFrame containing the data
                    x='stim_rec_dist',    # Column name for x-axis values
                    y=metric,             # Column name for y-axis values
                    style='ch_num_neural',  # Variable that determines marker style
                    color='black',       # Color for all points (overridden by hue)
                    size='dataset_stim_polarity',        # Variable that determines marker size
                    sizes=(150,50),       # Range of marker sizes (min, max)
                    markers=True,         # Enable different marker styles
                    alpha=0.7,           # Transparency (0-1)
                    s=75,                # Marker size
                    ax=ax                # Matplotlib axis to plot on
                )

                # Set larger font sizes
                ax.tick_params(labelsize=12)
                
                # Only show y-axis labels and ticks for leftmost plots
                if col == 0:
                    ax.set_ylabel(metric_labels[metric], fontsize=14)
                else:
                    ax.set_ylabel('')
                    ax.set_yticklabels([])
                
                ax.set_xlabel('Distance (mm)', fontsize=14)
                
                # Remove legend from all plots
                ax.legend().remove()
                
                # Set both x and y-axis limits
                ax.set_ylim(row_ylims[row])
                ax.set_xlim(global_xlims)
                
                ax.set_title(f'Date: {date}', fontsize=12)
                ax.grid(True)
        
        # Remove empty subplots
        for idx in range(len(dates), nrows * ncols):
            row = idx // ncols
            col = idx % ncols
            fig.delaxes(axes[row, col])
        
        # Add legend to the right of the last plot in the first row
        first_row_last_plot = min(ncols - 1, len(dates) - 1)
        # Use the last plot of the first row arbitrarily as the source of the
        # legend handles and labels, so this assumes that all plots have the
        # same legend items.
        handles, labels = axes[0, first_row_last_plot].get_legend_handles_labels()
        leg = fig.legend(handles, labels, title='Stim Polarity', 
                        bbox_to_anchor=(0.82, 0.5), loc='center left',
                        fontsize=12, title_fontsize=12)
        
        # Add title with proper spacing
        plt.suptitle(f'{nerve} - {metric_labels[metric]}', y=0.98, fontsize=16)
        # plt.tight_layout()
        
        if save_fig:
            plt.savefig(f'{self.output_path}/{metric}_{nerve}_grid.png', bbox_inches='tight', dpi=300)
        
        # Position the window on the primary monitor
        try:
            mngr = plt.get_current_fig_manager()
            mngr.window.wm_geometry("+0+0")  # Position at top-left
            # Set a reasonable max size that will fit on most monitors
            mngr.window.wm_geometry("1200x800")  # Width x Height in pixels
        except:
            pass  # In case this fails on some systems
    
        # save_file_name = f'{metric}_{nerve}_grid.png'
        # plt.savefig(save_file_name)
        plt.show()

    def plot_combined_dates(self, metric: str, nerve: str, fiber_type: str, save_fig: bool = False) -> None:
        """
        Create a single plot combining all experiment dates for a given metric and nerve type,
        using seaborn for enhanced visualization.
        
        Args:
            metric (str): The metric to plot ('peak_amplitude', 'peak_time', or 'auc')
            nerve (str): The nerve type to plot
            fiber_type (str): myel or unmyel
            save_fig (bool): Whether to save the figure
        """
        # Filter data for the specified nerve
        nerve_data = self.data_df[(self.data_df['rec_nerve'] == nerve) & (self.data_df['target_fiber_types'] == fiber_type)].copy()
        
        # Commented this out because it is no longer true!
        # Print warning that only myelinated fibers are being plotted 
        # print(f"Warning: Only myelinated fibers are being plotted for {nerve}.")
        
        # Dictionary for metric labels
        metric_labels = {
            'peak_amplitude': 'Peak Amplitude (µV)',
            'peak_time': 'Peak Time (ms)',
            'auc': 'AUC'
        }
        
        # Create figure with appropriate size
        plt.figure(figsize=(12, 8))
        
        # Create the scatter plot using seaborn
        sns.scatterplot(
            data=nerve_data,
            x='stim_rec_dist',
            y=metric,
            hue='expdate',           # Color by experiment date
            style='ch_num_neural',   # Different markers for channels
            size='dataset_stim_polarity',  # Size by polarity
            sizes=(100, 200),        # Range of marker sizes
            alpha=0.7,
            palette='husl',          # Use a colorful palette for dates
        )
        
        # Customize the plot
        plt.xlabel('Distance (mm)', fontsize=12)
        plt.ylabel(metric_labels[metric], fontsize=12)
        plt.title(f'{nerve} - {metric_labels[metric]} (All Dates Combined) - {fiber_type}', fontsize=14)
        
        # Customize legend
        plt.legend(
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            borderaxespad=0,
            fontsize=10
        )
        
        # Add grid with light opacity
        plt.grid(True, alpha=0.3)
        
        # Adjust layout to prevent legend cutoff
        plt.tight_layout()
        
        if save_fig:
            plt.savefig(f'{self.output_path}/{metric}_{nerve}_combined_dates_seaborn.png', 
                    bbox_inches='tight', dpi=300)
        
        plt.show()

# --- Load configuration ---
CONFIG_FILE = "local_config.yml"

try:
    with open(CONFIG_FILE, 'r') as f:
        # Use yaml.safe_load() for security
        local_config = yaml.safe_load(f) 
except FileNotFoundError:
    raise FileNotFoundError(f"Config file '{CONFIG_FILE}' not found. "
                          f"Did you copy 'config.yml.example' to 'config.yml'?")

RAW_DATA_PATH = local_config.get("raw_data_path")
PREPROCCESSED_CAP_PATH = local_config.get("preprocessed_CAP_path")
LOCAL_OUTPUT_PATH = local_config.get("local_output_path")

# Assert that the paths exist
assert Path(RAW_DATA_PATH).exists(), f"Data path '{RAW_DATA_PATH}' does not exist. Please check your local_config.yml file."
assert Path(PREPROCCESSED_CAP_PATH).exists(), f"Preprocessed CAP path '{PREPROCCESSED_CAP_PATH}' does not exist. Please check your local_config.yml file."
assert Path(LOCAL_OUTPUT_PATH).exists(), f"Local output path '{LOCAL_OUTPUT_PATH}' does not exist. Please check your local_config.yml file."

# --- End of config ---
# load df into a DataProcessor obj (if it exists)
STORED_DF_PATH = os.path.join(LOCAL_OUTPUT_PATH,"stored_df.pkl")
if Path(STORED_DF_PATH).exists():
    stored_df = pd.read_pickle(STORED_DF_PATH)
    processor = CNAPDataProcessor(str(Path(RAW_DATA_PATH) / 'AnalysisScript_CAPs_Batch_parameters.xlsx'),LOCAL_OUTPUT_PATH, stored_df)
else:
    # Load the preprocessed .mat files into one big DataFrame, but add some extra
    # features from the CAPs (e.g., auc, peak; see
    # CNAPDataProcessor.load_mat_file()). There is one folder of .mat files per row
    # in the Excel metadata file. The resulting DataFrame is in processor.data_df
    # processor = CNAPDataProcessor('AnalysisScript_CAPs_Batch_parameters.xlsx')
    processor = CNAPDataProcessor(str(Path(RAW_DATA_PATH) / 'AnalysisScript_CAPs_Batch_parameters.xlsx'),LOCAL_OUTPUT_PATH)
    processor.process_data(PREPROCCESSED_CAP_PATH, flag_plot=False) #FlagPlot to plot the gutcheck bounds
    
    #Save DataFrame to a .pkl file (if it doesn't already exist)
    processor.data_df.to_pickle(STORED_DF_PATH)


# Plot summary plots: plot_metric_grid plots all the features below from each
# nerve as a separate subplot; plot_combined_dates combines all the features
# from all nerves into a single plots
plot_summary_myel = False
if plot_summary_myel:
    fiber_type = 'myel'
    nerve_str_tuple = ('cVN','SN','tibial')
    for nerve_str in nerve_str_tuple:
        for metric in ['peak_amplitude', 'peak_time', 'auc']:
            processor.plot_combined_dates(metric, nerve=nerve_str, fiber_type=fiber_type)
            processor.plot_metric_grid(metric, nerve=nerve_str, fiber_type=fiber_type, ncols=8)

# Plot all the CAPs for each unique exp date that is from cervical vagus nerve
plot_CAPS_myel = True
if plot_CAPS_myel:
    nerve_str_tuple = ('cVN','SN','tibial')
    fiber_type = 'myel'
    for nerve_str in nerve_str_tuple:
        matching_expdates = sorted(processor.data_df[(processor.data_df['rec_nerve'] == nerve_str) & (processor.data_df['target_fiber_types'] == fiber_type)]['expdate'].unique())
        for exp_date in matching_expdates:
            processor.plot_data(exp_date=exp_date, fiber_type=fiber_type, nerve=nerve_str, show_analysis=True)

# TODO: The stim artifact window is not set correctly for unmyelinated fibers yet.
# Same plots as the above two sections but for unmyel
plot_summary_unmyel = False
if plot_summary_unmyel:
    fiber_type = 'unmyel'
    nerve_str = 'cVN'  # Only the cVN has unmyelinated fiber data
    for metric in ['peak_amplitude', 'peak_time', 'auc']:
        processor.plot_combined_dates(metric, nerve=nerve_str, fiber_type=fiber_type)
        processor.plot_metric_grid(metric, nerve=nerve_str, fiber_type=fiber_type, ncols=8)

plot_CAPS_unmyel = False
if plot_CAPS_unmyel:
    fiber_type = 'unmyel'
    nerve_str = 'cVN' # Only the cVN has unmyelinated fiber data
    matching_expdates = sorted(processor.data_df[(processor.data_df['rec_nerve'] == nerve_str) & (processor.data_df['target_fiber_types'] == fiber_type)]['expdate'].unique())
    for exp_date in matching_expdates:
        processor.plot_data(exp_date=exp_date, fiber_type=fiber_type, nerve=nerve_str, show_analysis=True)
    

