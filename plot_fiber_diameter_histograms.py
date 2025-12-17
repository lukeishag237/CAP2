import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import yaml
import os

def get_color_mapping(df: pd.DataFrame) -> dict:
    """
    Create a consistent color mapping for all nerve names in the dataset,
    using Set1, Set2, and Dark2 palettes for different nerve types.
    
    Args:
        df: DataFrame containing the nerve data
        
    Returns:
        dict: Mapping of nerve names to colors
    """
    # Get unique nerve types
    nerve_types = sorted(df['NerveType'].unique())
    
    # Specify palettes for each nerve type
    palette_map = {
        nerve_types[0]: 'Dark2',
        nerve_types[1]: 'Set2',
        nerve_types[2]: 'Set1'
    }
    
    color_dict = {}
    
    # Create color mapping for each nerve type
    for nerve_type in nerve_types:
        # Get names for this nerve type
        names_this_type = sorted(df[df['NerveType'] == nerve_type]['Name'].unique())
        
        # Get colors for this nerve type
        colors = sns.color_palette(palette_map[nerve_type], n_colors=len(names_this_type))
        
        # Add to color dictionary
        color_dict.update(dict(zip(names_this_type, colors)))
    
    return color_dict

def plot_diameter_histograms_grid(csv_file: str, nerve_type: str, ncols: int = 3, 
                                save_fig: bool = False) -> None:
    """
    Create a grid of histograms showing fiber diameter distributions for each nerve
    of a specified nerve type, using seaborn's FacetGrid.
    
    Args:
        csv_file (str): Path to the CSV file containing the data
        nerve_type (str): Type of nerve to analyze
        ncols (int): Number of columns in the grid
        save_fig (bool): Whether to save the figure
    """
    # Load the data
    df = pd.read_csv(csv_file)
    
    # Get consistent color mapping
    color_dict = get_color_mapping(df)
    
    # Filter for the specified nerve type
    nerves_df = df[df['NerveType'] == nerve_type]
    
    # Calculate bin edges with 0.07 μm resolution
    # Find the maximum diameter value to set upper bound
    fiber_diameter_resolution = 0.07  # μm ; based on Fig 9 of (Peña et al., 2024, PMID: 38427699)
    max_diameter = nerves_df['MinFeret'].max()
    bin_edges = np.arange(0, max_diameter + 0.07, 0.07)  # From 0 to max + 0.07 in 0.07 μm steps
    
    # Calculate number of rows needed
    n_names = len(nerves_df['Name'].unique())
    nrows = (n_names + ncols - 1) // ncols
    
    # Create FacetGrid
    g = sns.FacetGrid(
        nerves_df,
        col='Name',
        col_wrap=ncols,
        height=5,
        aspect=0.8,
        sharey=True  # Share y axis scales
    )
    
    # Add histograms to each facet
    g.map_dataframe(
        sns.histplot,
        x='MinFeret',
        bins=bin_edges,  # Use explicit bin edges for 0.07 μm resolution
        kde=True,
        kde_kws={
            'bw_adjust': 0.5,  # Adjust KDE bandwidth; 1=default, <1=rougher, >1=smoother
        }
    )
    
    # Customize using custom colors
    for ax, name in zip(g.axes.flat, sorted(nerves_df['Name'].unique())):
        for patch in ax.patches:
            patch.set_facecolor('lightgray')  # Make histogram bars light gray for visibility
            patch.set_edgecolor('none')  # Remove edge colors
        for line in ax.lines:
            line.set_color(color_dict[name])  
    
    # Add super title
    plt.suptitle(f'Fiber Diameter Distributions - {nerve_type}', 
                y=0.95, fontsize=12)
    
    # Set axis labels
    g.set_axis_labels('Fiber Diameter (μm)', 'Count')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure; force overwrite if it exists
    if save_fig:
        plt.savefig(f'diameter_histograms_{nerve_type}.svg', 
                   bbox_inches='tight', dpi=300)
        
    
    plt.show()


def plot_diameter_histograms_grid_with_gmm(csv_file: str, nerve_type: str, ncols: int = 3, 
                                         save_fig: bool = False) -> dict:
    """
    Create a grid of histograms showing fiber diameter distributions for each nerve
    of a specified nerve type, with fitted 2-component Gaussian Mixture Models overlaid.
    
    Args:
        csv_file (str): Path to the CSV file containing the data
        nerve_type (str): Type of nerve to analyze
        ncols (int): Number of columns in the grid
        save_fig (bool): Whether to save the figure
        
    Returns:
        dict: GMM parameters for each nerve (means, std devs, weights, log-likelihood)
    """
    # Load the data
    df = pd.read_csv(csv_file)
    
    # Get consistent color mapping
    color_dict = get_color_mapping(df)
    
    # Filter for the specified nerve type
    nerves_df = df[df['NerveType'] == nerve_type]
    
    # Calculate bin edges with 0.07 μm resolution
    fiber_diameter_resolution = 0.07  # μm
    max_diameter = nerves_df['MinFeret'].max()
    bin_edges = np.arange(0, max_diameter + 0.07, 0.07)
    
    # Calculate number of rows needed
    nerve_names = sorted(nerves_df['Name'].unique())
    n_names = len(nerve_names)
    nrows = (n_names + ncols - 1) // ncols
    
    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*5, nrows*5), 
                            sharey=True, squeeze=False)
    axes = axes.flatten()
    
    # Dictionary to store GMM results
    gmm_results = {}
    
    # Plot each nerve
    for i, name in enumerate(nerve_names):
        ax = axes[i]
        
        # Get data for this nerve
        nerve_data = nerves_df[nerves_df['Name'] == name]['MinFeret'].values
        
        # Create histogram
        counts, bin_edges_hist, patches = ax.hist(nerve_data, bins=bin_edges, 
                                                 alpha=0.7, color='lightgray', 
                                                 edgecolor='none',
                                                 density=False)

        # Fit 2-component GMM
        try:
            gmm = GaussianMixture(n_components=2, random_state=42)
            gmm.fit(nerve_data.reshape(-1, 1))
            
            # Extract parameters
            means = gmm.means_.flatten()
            covariances = gmm.covariances_.flatten()
            std_devs = np.sqrt(covariances)
            weights = gmm.weights_
            log_likelihood = gmm.score(nerve_data.reshape(-1, 1))
            
            # Store results
            gmm_results[name] = {
                'means': means,
                'std_devs': std_devs,
                'weights': weights,
                'log_likelihood': log_likelihood,
                'n_fibers': len(nerve_data)
            }
            
            # Create smooth x values for plotting
            x_smooth = np.linspace(0, max_diameter, 1000)

            # Calculate bin width
            bin_width = bin_edges[1] - bin_edges[0]

            # Calculate overall mixture density
            mixture_density = np.zeros_like(x_smooth)
            for j in range(2):
                component_density = weights[j] * norm.pdf(x_smooth, means[j], std_devs[j])
                mixture_density += component_density
                
                # Plot individual components (dashed lines)
                ax.plot(x_smooth, len(nerve_data) * component_density * bin_width, '--', 
                       color=color_dict[name], alpha=0.6, linewidth=1.5,
                       label=f'Component {j+1}' if i == 0 else "")
            
            # Plot overall mixture (solid line)
            ax.plot(x_smooth, len(nerve_data) * mixture_density * bin_width, '-', 
                   color=color_dict[name], linewidth=2.5,
                   label='GMM Fit' if i == 0 else "")
            
        except Exception as e:
            print(f"Warning: GMM fitting failed for {name}: {e}")
            gmm_results[name] = None
        
        # Customize subplot
        ax.set_title(f'{name}', fontsize=12)
        ax.set_xlabel('Fiber Diameter (μm)')
        if i % ncols == 0:  # First column
            ax.set_ylabel('Density')
        
        # Set consistent x-axis limits
        ax.set_xlim(0, 18)
    
    # Hide unused subplots
    for i in range(n_names, len(axes)):
        axes[i].set_visible(False)
    
    # Add overall title
    fig.suptitle(f'Fiber Diameter Distributions with GMM Fits - {nerve_type}', 
                fontsize=14, y=0.95)
    
    # Add legend to first subplot
    if len(nerve_names) > 0:
        axes[0].legend(loc='upper right', fontsize=10)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    if save_fig:
        plt.savefig(f'diameter_histograms_gmm_{nerve_type}.svg', 
                   bbox_inches='tight', dpi=300)
    
    plt.show()
    
    # Print GMM results summary
    print(f"\nGMM Fitting Results for {nerve_type}:")
    print("=" * 60)
    for name, results in gmm_results.items():
        if results is not None:
            print(f"\n{name} (n={results['n_fibers']} fibers):")
            print(f"  Component 1: μ={results['means'][0]:.2f}μm, σ={results['std_devs'][0]:.2f}μm, weight={results['weights'][0]:.3f}")
            print(f"  Component 2: μ={results['means'][1]:.2f}μm, σ={results['std_devs'][1]:.2f}μm, weight={results['weights'][1]:.3f}")
            print(f"  Log-likelihood: {results['log_likelihood']:.2f}")
        else:
            print(f"\n{name}: GMM fitting failed")

        # Print Markdown/CSV-style table for easy copy-paste with summary and legend
        print("\nCopy-paste table for summary (tab-delimited):")
        print("Nerve\tA\tB\tC\tD\tE\tF")

        nerves = []
        n_fibers = []
        gmm1_mean = []
        gmm1_cv = []
        gmm2_mean = []
        gmm2_cv = []
        gmm1_weight = []

        for name, results in gmm_results.items():
            if results is not None:
                cv1 = results['std_devs'][0] / results['means'][0] if results['means'][0] != 0 else float('nan')
                cv2 = results['std_devs'][1] / results['means'][1] if results['means'][1] != 0 else float('nan')
                print(f"{name}\t{results['n_fibers']}\t{results['means'][0]:.2f}\t{cv1:.2f}\t{results['means'][1]:.2f}\t{cv2:.2f}\t{results['weights'][0]:.3f}")
                nerves.append(name)
                n_fibers.append(results['n_fibers'])
                gmm1_mean.append(results['means'][0])
                gmm1_cv.append(cv1)
                gmm2_mean.append(results['means'][1])
                gmm2_cv.append(cv2)
                gmm1_weight.append(results['weights'][0])
            else:
                print(f"{name}\t---\t---\t---\t---\t---\t---")

        def mean_str(arr, fmt=".2f"):
            if len(arr) == 0:
                return "---"
            arr = np.array(arr)
            return f"{arr.mean():{fmt}}"
        def range_str(arr, fmt=".2f"):
            if len(arr) == 0:
                return "---"
            arr = np.array(arr)
            return f"[{arr.min():{fmt}}-{arr.max():{fmt}}]"

        print(f"Mean\t{mean_str(n_fibers, '.0f')}\t{mean_str(gmm1_mean)}\t{mean_str(gmm1_cv)}\t{mean_str(gmm2_mean)}\t{mean_str(gmm2_cv)}\t{mean_str(gmm1_weight, '.3f')}")
        print(f"Range\t{range_str(n_fibers, '.0f')}\t{range_str(gmm1_mean)}\t{range_str(gmm1_cv)}\t{range_str(gmm2_mean)}\t{range_str(gmm2_cv)}\t{range_str(gmm1_weight, '.3f')}")

        print("\nLegend:")
        print("A: Number of Fibers per nerve")
        print("B: Gaussian Mixture Model Component 1 (Mean, μm)")
        print("C: Gaussian Mixture Model Component 1 (Std/Mean, CV)")
        print("D: Gaussian Mixture Model Component 2 (Mean, μm)")
        print("E: Gaussian Mixture Model Component 2 (Std/Mean, CV)")
        print("F: Mix Proportion of Component 1")
    
    return gmm_results


def plot_combined_diameter_distributions(csv_file: str, save_fig: bool = False) -> None:
    """
    Create separate plots for each nerve type, combining all nerves of that type
    into a single plot with different colors.
    
    Args:
        csv_file (str): Path to the CSV file containing the data
        save_fig (bool): Whether to save the figure
    """
    # Load the data
    df = pd.read_csv(csv_file)
    
    # Get consistent color mapping
    color_dict = get_color_mapping(df)
    
    # Get unique nerve types
    nerve_types = sorted(df['NerveType'].unique())
    
    # Create a separate plot for each nerve type
    for nerve_type in nerve_types:
        # Filter data for this nerve type
        nerve_type_df = df[df['NerveType'] == nerve_type]
        
        # Set font size for better readability
        plt.rcParams['font.size'] = 14 

        # Create new figure
        plt.figure()
        
        # Plot all distributions at once using seaborn's built-in functionality
        sns.kdeplot(
            data=nerve_type_df,
            x='MinFeret',
            hue='Name',
            palette=color_dict,
            bw_adjust=0.5,  # Adjust bandwidth for smoother curves; 1 means default; less than 1 makes it rougher, more than 1 makes it smoother
            alpha=0.9,
            common_norm=False,  # Each nerve's KDE area equals its fiber count
        )

        # Customize plot
        plt.title(f'Fiber Diameter Distributions - {nerve_type}')
        plt.xlabel('Fiber Diameter (μm)')
        # Label y-axis appropriately since KDE area now equals fiber count
        plt.ylabel('Density (Fibers per μm)')
        
        plt.xlim(0, 18)  # Set x-axis limit to 18 um for better visualization

        if save_fig:
            plt.savefig(f'diameter_distributions_{nerve_type}.svg', 
                       bbox_inches='tight', dpi=300)
        
        plt.show()

def plot_fiber_counts(csv_file: str, save_fig: bool = False) -> pd.DataFrame:
    """
    Create a scatter plot of fiber counts by nerve type, with means and standard deviations.
    Points are colored by nerve name, with mean and std shown as error bars.
    
    Args:
        csv_file (str): Path to the CSV file containing the data
        save_fig (bool): Whether to save the figure
        
    Returns:
        pd.DataFrame: Summary statistics for each nerve type
    """
    # Load the data
    df = pd.read_csv(csv_file)
    
    # Get consistent color mapping
    color_dict = get_color_mapping(df)
    
    # Count fibers for each combination of NerveType and Name
    fiber_counts = df.groupby(['NerveType', 'Name']).size().reset_index(name='Fiber_Count')
    
    # Calculate summary statistics for each nerve type
    summary_stats = fiber_counts.groupby('NerveType').agg({
        'Fiber_Count': ['mean', 'std', 'count']
    }).round(2)
    summary_stats.columns = ['Mean_Fibers', 'Std_Fibers', 'N_Samples']
    
    fiber_diameters_mean = df.groupby(['NerveType', 'Name'])['MinFeret'].mean().reset_index(name='Mean_Diameter')
    fiber_diameters_std = df.groupby(['NerveType', 'Name'])['MinFeret'].std().reset_index(name='Std_Diameter')
    fiber_diameters_cv = fiber_diameters_std.copy()
    fiber_diameters_cv['CV_Diameter'] = fiber_diameters_std['Std_Diameter'] / fiber_diameters_mean['Mean_Diameter']
    fiber_diameters_range = df.groupby(['NerveType', 'Name'])['MinFeret'].agg(['min', 'max']).reset_index()

    # Print fiber diameter stats for each nerve
    print("\nFiber Diameter Statistics by Nerve:")
    print("----------------------------------")

    for name in sorted(df['Name'].unique()):
        mean_dia = fiber_diameters_mean[fiber_diameters_mean['Name'] == name]['Mean_Diameter'].values[0]
        std_dia = fiber_diameters_std[fiber_diameters_std['Name'] == name]['Std_Diameter'].values[0]
        cv_dia = fiber_diameters_cv[fiber_diameters_cv['Name'] == name]['CV_Diameter'].values[0]
        dia_min = fiber_diameters_range[fiber_diameters_range['Name'] == name]['min'].values[0]
        dia_max = fiber_diameters_range[fiber_diameters_range['Name'] == name]['max'].values[0]
        print(f"{name}: Mean Diameter = {mean_dia:.2f} μm, Std Dev = {std_dia:.2f} μm, CV = {cv_dia:.2f}, Range = [{dia_min:.2f}-{dia_max:.2f}] μm")

    # Print fiber diameter stats for each nerve as above, but grouped by nerve type
    print("\nFiber Diameter Statistics by Nerve Type:")
    print("---------------------------------------")
    for nerve_type in summary_stats.index:
        nerves_of_type = sorted(df[df['NerveType'] == nerve_type]['Name'].unique())
        print(f"\n{nerve_type}:")
        # Print header for tab-delimited format
        print("  Nerve\tFiber Count\tMean (μm)\tStd/Mean\tRange (μm)")
        mean_arrays = np.array([])
        std_arrays = np.array([])
        cv_arrays = np.array([])
        fiber_count_arrays = np.array([])
        for name in nerves_of_type:
            mean_dia = fiber_diameters_mean[fiber_diameters_mean['Name'] == name]['Mean_Diameter'].values[0]
            std_dia = fiber_diameters_std[fiber_diameters_std['Name'] == name]['Std_Diameter'].values[0]
            cv_dia = fiber_diameters_cv[fiber_diameters_cv['Name'] == name]['CV_Diameter'].values[0]
            dia_min = fiber_diameters_range[fiber_diameters_range['Name'] == name]['min'].values[0]
            dia_max = fiber_diameters_range[fiber_diameters_range['Name'] == name]['max'].values[0]
            count = fiber_counts[fiber_counts['Name'] == name]['Fiber_Count'].values[0]
            # print(f"  {name}: Mean Diameter = {mean_dia:.2f} μm, Std Dev = {std_dia:.2f} μm, CV = {cv_dia:.2f}, Range = [{dia_min:.2f}-{dia_max:.2f}] μm")
            # as above, but use tab-delimited format for easy copy-paste
            print(f"  {name}\t{count}\t{mean_dia:.2f}\t{cv_dia:.2f}\t[{dia_min:.2f}-{dia_max:.2f}]")
            mean_arrays = np.append(mean_arrays, mean_dia)
            std_arrays = np.append(std_arrays, std_dia)
            cv_arrays = np.append(cv_arrays, cv_dia)
            fiber_count_arrays = np.append(fiber_count_arrays, count)

        # Print mean for each nerve type in tab delimited format as a single row (e.g., average of mean_arrays, cv_arrays, and fiber_count_arrays)
        print(f" Mean\t{np.mean(fiber_count_arrays):.0f}\t{np.mean(mean_arrays):.2f}\t{np.mean(cv_arrays):.2f}")

        # Print the range for each nerve type in tab delimited format as a single row (e.g., min and max of mean_arrays, cv_arrays, and fiber_count_arrays)
        print(f" Range\t[{np.min(fiber_count_arrays):.0f}-{np.max(fiber_count_arrays):.0f}]\t[{np.min(mean_arrays):.2f}-{np.max(mean_arrays):.2f}]\t[{np.min(cv_arrays):.2f}-{np.max(cv_arrays):.2f}]")

    plt.rcParams['font.size'] = 14 

    # Create the plot
    plt.figure()
    
    # Create scatter plot using swarmplot with custom palette
    sns.swarmplot(
        data=fiber_counts,
        x='NerveType',
        y='Fiber_Count',
        hue='Name',
        size=8,
        alpha=1,
        palette=color_dict
    )
    
    # Add mean and standard deviation
    for i, nerve_type in enumerate(summary_stats.index):
        stats = summary_stats.loc[nerve_type]
        mean = stats['Mean_Fibers']
        std = stats['Std_Fibers']
        
        # Plot mean as a horizontal line
        plt.hlines(y=mean, xmin=i-0.1, xmax=i+0.1, color='black', linewidth=1, label='Mean' if i==0 else "")
        
        # Plot standard deviation as vertical error bars
        plt.errorbar(x=i, y=mean, yerr=std, color='black', capsize=5, 
                    capthick=1, linewidth=1, label='Std Dev' if i==0 else "")
    
    # Customize plot
    plt.title('Fiber Counts by Nerve Type')
    plt.xlabel('Nerve Type')
    plt.ylabel('Number of Fibers')
    
    # Rotate x-axis labels if needed
    plt.xticks(rotation=45)
    
    # Add legend with both scatter points and statistics
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(
        handles=handles,
        labels=labels,
        title='Nerve Name',
        bbox_to_anchor=(1.05, 1),
        loc='upper left'
    )
    
    # Adjust layout
    plt.tight_layout()
    
    if save_fig:
        plt.savefig('fiber_counts_scatter.svg', bbox_inches='tight', dpi=300)
    
    plt.show()
    
    # Print the fiber counts for each nerve and the mean fiber diameter for each nerve and the std/mean (CV) and range of fiber diameters
    print("\nFiber Counts by Nerve:")
    print("-----------------------")
    for _, row in fiber_counts.iterrows():
        print(f"{row['Name']} ({row['NerveType']}): {row['Fiber_Count']} fibers")


    # Print summary statistics
    print("\nSummary Statistics by Nerve Type:")
    print("----------------------------------")
    for nerve_type in summary_stats.index:
        stats = summary_stats.loc[nerve_type]
        print(f"\n{nerve_type}:")
        print(f"  Mean Fiber Count: {stats['Mean_Fibers']:.2f}")
        print(f"  Standard Deviation: {stats['Std_Fibers']:.2f}")
        print(f"  Number of Samples: {stats['N_Samples']}")
    
    return summary_stats


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

csv_file = os.path.join(RAW_DATA_PATH, '..', '60x', 'derivatives', 'summary', 'all_samples_fiber_diameters.csv')

stats = plot_fiber_counts(csv_file, save_fig=True)
print(stats)

plot_diameter_histograms_grid(csv_file, 'sciatic', save_fig=True, ncols=5)

plot_diameter_histograms_grid(csv_file, 'tibial', save_fig=True, ncols=5)

plot_diameter_histograms_grid(csv_file, 'vagus', save_fig=True, ncols=5)

# Test the new GMM-based function
gmm_results = plot_diameter_histograms_grid_with_gmm(csv_file, 'vagus', save_fig=True, ncols=5)

plot_combined_diameter_distributions(csv_file, save_fig=True)   

# plot_diameter_histograms_grid('C:\\Users\\hroth\\Box\\adhoc_grill_lab_histology\\Edgar Pena\\60x\\derivatives\\fiber_stats\\auto_contour_stats_with_SAE.csv', 'sciatic', save_fig=True, ncols=5)
# plot_diameter_histograms_grid('C:\\Users\\hroth\\Box\\adhoc_grill_lab_histology\\Edgar Pena\\60x\\derivatives\\fiber_stats\\auto_contour_stats_with_SAE.csv', 'tibial', save_fig=True, ncols=5)
