import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from pathlib import Path
import numpy as np

def read_bed_file(filepath):
    """
    Read a BED file and return a DataFrame with peak information.
    Assumes 3-column BED format: chr, start, end
    """
    try:
        # Read BED file - only 3 columns (chr, start, end)
        df = pd.read_csv(filepath, sep='\t', header=None, 
                        names=['chr', 'start', 'end'])
        
        # Ensure start and end are numeric
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        
        # Remove any rows with non-numeric coordinates
        df = df.dropna(subset=['start', 'end'])
        
        # Calculate peak length
        df['length'] = df['end'] - df['start']
        
        return df
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return pd.DataFrame()

def plot_combined_peak_distribution(bed_directory, output_file=None):
    """
    Plot a single combined distribution of ALL peaks from ALL BED files.
    """
    # Find all BED files in directory
    bed_files = glob.glob(os.path.join(bed_directory, "*.bed"))
    
    if not bed_files:
        print(f"No BED files found in {bed_directory}")
        return
    
    print(f"Found {len(bed_files)} BED files")
    
    # Read all BED files and combine
    all_peaks = []
    total_peaks = 0
    
    for bed_file in bed_files:
        df = read_bed_file(bed_file)
        if not df.empty:
            all_peaks.append(df)
            print(f"Loaded {len(df)} peaks from {os.path.basename(bed_file)}")
            total_peaks += len(df)
    
    if not all_peaks:
        print("No valid BED files could be read")
        return
    
    # Combine ALL peaks into one dataset
    combined_df = pd.concat(all_peaks, ignore_index=True)
    
    # Filter out invalid peaks (negative or zero length)
    combined_df = combined_df[combined_df['length'] > 0]
    
    print(f"Total valid peaks from all files: {len(combined_df):,}")
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Single histogram of ALL peak lengths
    plt.hist(combined_df['length'], bins=50, alpha=0.7, edgecolor='black', color='skyblue')
    
    plt.xlabel('Peak Length (bp)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title(f'Combined Peak Length Distribution\n({len(combined_df):,} peaks from {len(bed_files)} files)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add some statistics as text on the plot
    mean_length = combined_df['length'].mean()
    median_length = combined_df['length'].median()
    
    stats_text = f'Mean: {mean_length:.1f} bp\nMedian: {median_length:.1f} bp\nMin: {combined_df["length"].min()} bp\nMax: {combined_df["length"].max()} bp'
    plt.text(0.7, 0.8, stats_text, transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
             verticalalignment='top', fontsize=10)
    
    plt.tight_layout()
    
    # Save plot if output file specified
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*50)
    print("SUMMARY STATISTICS")
    print("="*50)
    print(f"Total number of peaks: {len(combined_df):,}")
    print(f"Number of files processed: {len(bed_files)}")
    print(f"Average peak length: {combined_df['length'].mean():.1f} bp")
    print(f"Median peak length: {combined_df['length'].median():.1f} bp")
    print(f"Standard deviation: {combined_df['length'].std():.1f} bp")
    print(f"Peak length range: {combined_df['length'].min()} - {combined_df['length'].max()} bp")
    print(f"Number of chromosomes: {combined_df['chr'].nunique()}")
    
    return combined_df

# Alternative plotting styles
def plot_with_kde(bed_directory, output_file=None):
    """
    Create histogram + KDE overlay for combined peak distribution
    """
    bed_files = glob.glob(os.path.join(bed_directory, "*.bed"))
    all_peaks = []
    
    for bed_file in bed_files:
        df = read_bed_file(bed_file)
        if not df.empty:
            all_peaks.append(df)
    
    combined_df = pd.concat(all_peaks, ignore_index=True)
    combined_df = combined_df[combined_df['length'] > 0]
    
    plt.figure(figsize=(12, 8))
    
    # Histogram + KDE
    plt.hist(combined_df['length'], bins=50, alpha=0.6, density=True, 
             color='lightblue', edgecolor='black', label='Histogram')
    
    # Add KDE if there's variance
    if combined_df['length'].var() > 0:
        sns.kdeplot(combined_df['length'], color='red', linewidth=2, label='KDE')
    
    plt.xlabel('Peak Length (bp)')
    plt.ylabel('Density')
    plt.title(f'Peak Length Distribution with KDE\n({len(combined_df):,} peaks)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    plt.show()

# Example usage
if __name__ == "__main__":
    # Specify your BED files directory
    bed_directory = "/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/new-merged-filtered"  # Change this to your directory
    
    # Optional: specify output file for saving the plot
    output_file = "combined_peak_distribution.png"
    
    # Generate the combined plot
    peak_data = plot_combined_peak_distribution(bed_directory, output_file)
    
    # Optional alternative with KDE overlay
    # plot_with_kde(bed_directory, "peak_distribution_kde.png")
    
    # Optional: Save combined data to CSV
    if peak_data is not None and not peak_data.empty:
        peak_data.to_csv("all_peaks_combined.csv", index=False)
        print("Combined peak data saved to all_peaks_combined.csv")