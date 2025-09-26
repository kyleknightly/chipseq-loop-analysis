#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, optimal_leaf_ordering

# Define custom colormap

def load_proportions(path):
    """Load TF proportions and total count"""
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    df["protein"] = df["protein"].str.strip()
    tf_to_prop = dict(zip(df["protein"], df["proportion"]))
    any_tf = df.iloc[0]
    total = any_tf["count"] / any_tf["proportion"]
    return tf_to_prop, total

def load_pair_counts(path):
    """Load TF pair counts from bash output"""
    df = pd.read_csv(path, sep="\t")
    return df

def build_matrix(pair_counts_df, all_tfs):
    """Build symmetric matrix from pair counts"""
    tf_to_idx = {tf: i for i, tf in enumerate(all_tfs)}
    n = len(all_tfs)
    matrix = np.zeros((n, n), dtype=int)
    
    for _, row in pair_counts_df.iterrows():
        tf1, tf2, count = row['tf1'], row['tf2'], row['count']
        if tf1 in tf_to_idx and tf2 in tf_to_idx:
            i, j = tf_to_idx[tf1], tf_to_idx[tf2]
            matrix[i, j] = count
            matrix[j, i] = count  # Symmetric
    
    return matrix

def plot_enrichment_heatmap(name, enrichment_df, distance_file=None):
    """
    Create clustered heatmap with dendrogram and optional distance-based label coloring
    """
    print(f"Creating plot: {name}")
    
    # Handle zeros by replacing with small values
    maxval = enrichment_df.max().max()
    log10_enrichment_df = np.log10(enrichment_df.replace(0, 10 ** (-np.log10(maxval))))

    vmin = log10_enrichment_df.min().min()
    vmax = log10_enrichment_df.max().max()

    # Perform hierarchical clustering
    linkage_matrix = linkage(enrichment_df, method='centroid')
    linkage_matrix_olo = optimal_leaf_ordering(linkage_matrix, enrichment_df)
    ordered_index = leaves_list(linkage_matrix_olo)
    ordered_protein_names = log10_enrichment_df.index[ordered_index].tolist()
    
    print(f"Ordered proteins for {name}: {ordered_protein_names}") 
    
    # Reorder the matrix
    ordered_log10_enrichment_df = log10_enrichment_df.iloc[ordered_index, ordered_index]

    # Create figure with dendrogram and heatmap
    fig, (ax_dendro, ax_heatmap) = plt.subplots(1, 2, figsize=(70, 50), 
                                               gridspec_kw={'width_ratios': [1, 15]})
    ax_heatmap.set_aspect('equal')
    
    # Plot dendrogram
    dendro = dendrogram(linkage_matrix_olo, labels=ordered_log10_enrichment_df.index, 
                       orientation='left', ax=ax_dendro)
    ax_dendro.invert_yaxis()
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    # Plot heatmap
    heatmap = sns.heatmap(
        ordered_log10_enrichment_df, 
        cmap='RdBu_r',
        annot=False,      
        linewidths=0.05,   
        center=0,  
        ax=ax_heatmap,
        xticklabels=True,
        yticklabels=True,
        vmin=vmin,  
        vmax=vmax   
    )

    # Color labels based on distance to TSS if provided
    if distance_file and sys.path[0]:  # Check if file exists
        try:
            distance_df = pd.read_csv(distance_file,sep='\t', names =['protein', 'score'])
            distance_dict = dict(zip(distance_df["protein"], distance_df["score"]))
            
            # Get distances for normalization
            distances = [distance_dict[prot] for prot in enrichment_df.index 
                        if prot in distance_dict and not pd.isna(distance_dict[prot])]
            
            if distances:  # If we have distance data
                vmin_dist, vmax_dist = min(distances), max(distances)
                norm = mcolors.Normalize(vmin=vmin_dist, vmax=vmax_dist)
                cmap = mpl.colormaps["viridis"]
                
                # Color the labels
                for label in ax_heatmap.get_yticklabels():
                    protein = label.get_text()
                    distance = distance_dict.get(protein)
                    
                    if distance is not None and not pd.isna(distance):
                        color = cmap(norm(distance))
                        label.set_color(color)
                    else:
                        label.set_color("black")
                    label.set_fontsize(20)
            else:
                print("No valid distance data found")
        except Exception as e:
            print(f"Could not load distance file: {e}")
    
    # Stagger y-tick labels for better readability
    # for i, label in enumerate(ax_heatmap.get_yticklabels()):
    #     if i % 2 == 0:
    #         label.set_x(-0.0005)
    #     else:
    #         label.set_x(-0.018)

    # Adjust tick parameters
    ax_heatmap.yaxis.set_tick_params(which='both', length=0)
    # for i, tick in enumerate(ax_heatmap.yaxis.get_major_ticks()):
    #     tick.tick1line.set_markersize(2.75 if i % 2 == 0 else 47)

    # Set tick positions and labels
    heatmap.set_xticks(np.arange(len(ordered_log10_enrichment_df.columns)) + 0.5)
    heatmap.set_yticks(np.arange(len(ordered_log10_enrichment_df.index)) + 0.5)
    heatmap.set_xticklabels(ordered_log10_enrichment_df.columns, rotation=90, fontsize=9)
    heatmap.set_yticklabels(ordered_log10_enrichment_df.index, rotation=0, fontsize=9)

    # Format colorbar
    colorbar = heatmap.collections[0].colorbar
    colorbar.ax.yaxis.label.set_size(30)
    colorbar.ax.tick_params(labelsize=30)
    
    # Convert log10 ticks back to original scale
    # log_ticks = colorbar.get_ticks()
    # true_ticks = [10 ** tick for tick in log_ticks]
    # formatted_ticks = [f"{t:.2e}" if t < 1 else f"{int(t):,}" for t in true_ticks]
    # colorbar.set_ticks(log_ticks)
    # colorbar.set_ticklabels(formatted_ticks)

    plt.subplots_adjust(wspace=0.07)
    plt.savefig(name, dpi=300, bbox_inches='tight')
    plt.close()  # Free memory
    
    return ordered_protein_names

def main():
    EXCLUDE_PROTEINS = [
    'H3K4me1', 'H3K4me2', 'H3K4me3',   # Histone marks
    'H3K27ac', 'H3K27me3', 'H3K36me3',
    'H3K9ac', 'H3K9me3', 'H4K20me1',
    'H3K79me2', 'POLR2AphosphoS5'
    ]
    # File paths
    head = '/mnt/altnas/work/Kyle.Knightly/looppi/old-looppi/rand-1/out/'
    anchor_props_file = head+'04_tf_props.unique.tsv'
    end_props_file = head+'05_tf_props.weightedEnds.tsv'
    cis_counts_file = head+'06_cis_pairs.long.tsv'
    trans_counts_file = head+'08_trans_pairs.long.tsv'
    
    print("Loading proportions...")
    anchor_props, N_cis = load_proportions(anchor_props_file)
    print(N_cis) 
    end_props, N_trans = load_proportions(end_props_file)
    print(N_trans)
    
    # Get all TFs from proportion files
    all_tfs_raw = sorted(set(anchor_props.keys()) | set(end_props.keys()))
    all_tfs = [tf for tf in all_tfs_raw if tf not in EXCLUDE_PROTEINS]
    print(f"Found {len(all_tfs)} unique TFs")
    
    print("Loading pair counts...")
    cis_counts_df = load_pair_counts(cis_counts_file)
    trans_counts_df = load_pair_counts(trans_counts_file)
    
    print("Building matrices...")
    cis_obs = build_matrix(cis_counts_df, all_tfs)
    trans_obs = build_matrix(trans_counts_df, all_tfs)
    
    print("Computing expected matrices...")
    # Vectorized proportion arrays
    anchor_props_vec = np.array([anchor_props.get(tf, 0) for tf in all_tfs])
    end_props_vec = np.array([end_props.get(tf, 0) for tf in all_tfs])
    
    cis_exp = np.outer(anchor_props_vec, anchor_props_vec) * N_cis
    trans_exp = np.outer(end_props_vec, end_props_vec) * N_trans
    np.fill_diagonal(trans_exp, np.diag(trans_exp) / 2)
    
    print("Computing enrichment...")
    # Add pseudocounts to match pipeline.py: +1 to observed counts, +1 to expected counts
    cis_obs_pseudo = cis_obs + 1
    trans_obs_pseudo = trans_obs + 1
    cis_exp_pseudo = cis_exp + 1
    trans_exp_pseudo = trans_exp + 1
    
    cis_enrich = np.divide(cis_obs_pseudo, cis_exp_pseudo, out=np.zeros_like(cis_exp_pseudo, dtype=float), where=cis_exp_pseudo != 0)
    trans_enrich = np.divide(trans_obs_pseudo, trans_exp_pseudo, out=np.zeros_like(trans_exp_pseudo, dtype=float), where=trans_exp_pseudo != 0)
    
    print("Saving results...")
    pd.DataFrame(cis_enrich, index=all_tfs, columns=all_tfs).to_csv(head+'cis_enrichment_matrix.tsv', sep='\t')
    pd.DataFrame(trans_enrich, index=all_tfs, columns=all_tfs).to_csv(head+'trans_enrichment_matrix.tsv', sep='\t')
    
    # Filter for TFs with end_proportion >= 0.1
    print("Filtering for TFs with end_proportion >= 0.1...")
    filtered_tfs = [tf for tf in all_tfs if end_props.get(tf, 0) >= 0.03]
    print(f"Filtered from {len(all_tfs)} to {len(filtered_tfs)} TFs")
    
    # Get indices of filtered TFs in original matrix
    tf_to_idx = {tf: i for i, tf in enumerate(all_tfs)}
    filtered_indices = [tf_to_idx[tf] for tf in filtered_tfs]
    
    # Extract submatrices
    cis_filtered = cis_enrich[np.ix_(filtered_indices, filtered_indices)]
    trans_filtered = trans_enrich[np.ix_(filtered_indices, filtered_indices)]
    
    # Save filtered matrices
    pd.DataFrame(cis_filtered, index=filtered_tfs, columns=filtered_tfs).to_csv(head+'geq10p_cis_enrichment_matrix.tsv', sep='\t')
    pd.DataFrame(trans_filtered, index=filtered_tfs, columns=filtered_tfs).to_csv(head+'geq10p_trans_enrichment_matrix.tsv', sep='\t')
    
    print(f"✅ Done! Saved both full and filtered (>= 10%) enrichment matrices.")
    print(f"Full matrices shape: {cis_obs.shape}")
    print(f"Filtered matrices shape: {cis_filtered.shape}")
    
    
    print(f"Total CIS interactions: {cis_obs.sum()}")
    print(f"Total TRANS interactions: {trans_obs.sum()}")
    
    # Create plots
    print("Generating plots...")
    distance_file = "/mnt/altnas/work/Kyle.Knightly/tf_average_loop_sizes.tsv"
    
    # Plot full matrices
    plot_enrichment_heatmap(head+'cis_enrichment.png', 
                           pd.DataFrame(cis_enrich, index=all_tfs, columns=all_tfs),
                           distance_file)
    plot_enrichment_heatmap(head+'trans_enrichment.png', 
                           pd.DataFrame(trans_enrich, index=all_tfs, columns=all_tfs),
                           distance_file)

    plot_enrichment_heatmap(head+'geq10p_cis_enrichment.png',
                            pd.DataFrame(cis_filtered, index=filtered_tfs, columns=filtered_tfs),
                            distance_file)
    plot_enrichment_heatmap(head+'geq10p_trans_enrichment.png',
                            pd.DataFrame(trans_filtered, index=filtered_tfs, columns=filtered_tfs),
                            distance_file)
    
    print("🎨 All plots saved!")

if __name__ == "__main__":
    main()