#!/usr/bin/env python3
"""
Analyze network structure to diagnose motif analysis issues.
Check for dense connected components that might bias motif detection.
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def analyze_network_structure(edgelist_file):
    """Comprehensive analysis of network structure"""
    
    # Load network
    print("Loading network...")
    
    # Load network
    print("Loading network...")
    
    # Check file extension to determine format
    if edgelist_file.endswith('.gpickle') or edgelist_file.endswith('.pickle'):
        print("Loading pickled graph...")
        G = nx.read_gpickle(edgelist_file)
        print(f"Loaded pickled graph: {type(G)}")
    else:
        # First, peek at the file to understand format
        with open(edgelist_file, 'r') as f:
            first_line = f.readline().strip()
            columns = len(first_line.split())
        
        print(f"Detected {columns} columns in file")
        
        try:
            if columns == 2:
                # Two columns: source target
                G = nx.read_edgelist(edgelist_file, nodetype=int)
            elif columns == 3:
                # Three columns: source target weight
                G = nx.read_edgelist(edgelist_file, nodetype=int, data=(('weight', int),))
            else:
                # Try without specifying data format
                G = nx.read_edgelist(edgelist_file, nodetype=int, data=False)
        except Exception as e:
            print(f"Standard loading failed: {e}")
            print("Trying manual parsing...")
            
            # Manual parsing as fallback
            G = nx.Graph()  # or nx.DiGraph() if directed
            with open(edgelist_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            source, target = int(parts[0]), int(parts[1])
                            G.add_edge(source, target)
                        except ValueError:
                            print(f"Skipping line {line_num}: {line}")
                            continue
    
    print(f"Network loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Basic statistics
    print("\n=== BASIC NETWORK STATISTICS ===")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
    print(f"Density: {nx.density(G):.6f}")
    print(f"Average degree: {np.mean([d for n, d in G.degree()]):.2f}")
    
    # Connected components analysis
    print("\n=== CONNECTED COMPONENTS ANALYSIS ===")
    
    if G.is_directed():
        print("Directed graph detected - analyzing weakly connected components")
        components = list(nx.weakly_connected_components(G))
        strongly_connected = list(nx.strongly_connected_components(G))
        
        print(f"Weakly connected components: {len(components)}")
        print(f"Strongly connected components: {len(strongly_connected)}")
        
        # Use weakly connected components for main analysis
        component_sizes = [len(c) for c in components]
        strongly_component_sizes = [len(c) for c in strongly_connected]
        
        print(f"Largest weakly connected component: {max(component_sizes)}")
        print(f"Largest strongly connected component: {max(strongly_component_sizes)}")
        
    else:
        print("Undirected graph detected")
        components = list(nx.connected_components(G))
        component_sizes = [len(c) for c in components]
        
        print(f"Connected components: {len(components)}")
        print(f"Largest component size: {max(component_sizes)}")
    
    print(f"Smallest component size: {min(component_sizes)}")
    print(f"Average component size: {np.mean(component_sizes):.2f}")
    
    # Size distribution
    size_counts = Counter(component_sizes)
    print("\nComponent size distribution:")
    for size in sorted(size_counts.keys())[:10]:  # Show first 10 sizes
        count = size_counts[size]
        print(f"  Size {size}: {count} components")
    if len(size_counts) > 10:
        print(f"  ... ({len(size_counts)-10} more sizes)")
    
    # Analyze density within components
    print("\n=== COMPONENT DENSITY ANALYSIS ===")
    densities = []
    large_components = [c for c in components if len(c) >= 4]  # Only components that can have 4-node motifs
    
    for i, component in enumerate(large_components[:10]):  # Analyze first 10 large components
        subgraph = G.subgraph(component)
        density = nx.density(subgraph)
        densities.append(density)
        print(f"Component {i+1} (size {len(component)}): density = {density:.4f}")
    
    if densities:
        print(f"Average density within large components: {np.mean(densities):.4f}")
        print(f"Max density within components: {max(densities):.4f}")
        print(f"Min density within components: {min(densities):.4f}")
    
    # Clustering analysis
    print("\n=== CLUSTERING ANALYSIS ===")
    if G.is_directed():
        # For directed graphs, clustering is more complex
        try:
            clustering_coeffs = [nx.clustering(G, n) for n in G.nodes()]
            avg_clustering = np.mean(clustering_coeffs)
            print(f"Average clustering coefficient: {avg_clustering:.4f}")
            print(f"Max clustering coefficient: {max(clustering_coeffs):.4f}")
            
            # Compare to random directed network
            random_clustering = nx.average_clustering(nx.erdos_renyi_graph(G.number_of_nodes(), nx.density(G), directed=True))
            print(f"Random directed network clustering (same density): {random_clustering:.4f}")
            if random_clustering > 0:
                print(f"Clustering ratio (real/random): {avg_clustering/random_clustering:.2f}")
        except:
            print("Clustering calculation failed for directed graph")
    else:
        clustering_coeffs = [nx.clustering(G, n) for n in G.nodes()]
        avg_clustering = np.mean(clustering_coeffs)
        print(f"Average clustering coefficient: {avg_clustering:.4f}")
        print(f"Max clustering coefficient: {max(clustering_coeffs):.4f}")
        
        # Compare to random network
        random_clustering = nx.average_clustering(nx.erdos_renyi_graph(G.number_of_nodes(), nx.density(G)))
        print(f"Random network clustering (same density): {random_clustering:.4f}")
        print(f"Clustering ratio (real/random): {avg_clustering/random_clustering:.2f}")
    
    # Degree distribution analysis
    print("\n=== DEGREE DISTRIBUTION ===")
    if G.is_directed():
        in_degrees = [d for n, d in G.in_degree()]
        out_degrees = [d for n, d in G.out_degree()]
        total_degrees = [G.in_degree(n) + G.out_degree(n) for n in G.nodes()]
        
        print(f"In-degree statistics: min={min(in_degrees)}, max={max(in_degrees)}, mean={np.mean(in_degrees):.2f}")
        print(f"Out-degree statistics: min={min(out_degrees)}, max={max(out_degrees)}, mean={np.mean(out_degrees):.2f}")
        print(f"Total degree statistics: min={min(total_degrees)}, max={max(total_degrees)}, mean={np.mean(total_degrees):.2f}")
        
        # Check for hubs
        in_degree_threshold = np.mean(in_degrees) + 2 * np.std(in_degrees)
        out_degree_threshold = np.mean(out_degrees) + 2 * np.std(out_degrees)
        in_hubs = [n for n, d in G.in_degree() if d > in_degree_threshold]
        out_hubs = [n for n, d in G.out_degree() if d > out_degree_threshold]
        print(f"High in-degree nodes (>mean+2*std): {len(in_hubs)}")
        print(f"High out-degree nodes (>mean+2*std): {len(out_hubs)}")
        degrees = total_degrees
    else:
        degrees = [d for n, d in G.degree()]
        print(f"Degree statistics: min={min(degrees)}, max={max(degrees)}, mean={np.mean(degrees):.2f}, std={np.std(degrees):.2f}")
        
        # Check for hubs
        degree_threshold = np.mean(degrees) + 2 * np.std(degrees)
        hubs = [n for n, d in G.degree() if d > degree_threshold]
        print(f"Number of high-degree nodes (>mean+2*std): {len(hubs)}")
    
    # Motif analysis prediction
    print("\n=== MOTIF ANALYSIS PREDICTION ===")
    
    # Check if structure suggests motif bias
    bias_indicators = []
    
    if len(components) > G.number_of_nodes() * 0.1:  # Many small components
        bias_indicators.append("Many disconnected components")
    
    if avg_clustering > 0.3:  # High clustering
        bias_indicators.append("High clustering coefficient")
    
    if max(densities) if densities else 0 > nx.density(G) * 5:  # Some components much denser than average
        bias_indicators.append("Some components much denser than global network")
    
    if len(large_components) > 5 and np.std(component_sizes) > np.mean(component_sizes):
        bias_indicators.append("High variability in component sizes")
    
    if bias_indicators:
        print("⚠️  POTENTIAL MOTIF ANALYSIS BIAS DETECTED:")
        for indicator in bias_indicators:
            print(f"   - {indicator}")
        print("\nRecommendations:")
        print("   1. Analyze each large connected component separately")
        print("   2. Use degree-preserving randomization that maintains component structure")
        print("   3. Compare motif counts within components, not globally")
        print("   4. Consider using more sophisticated null models")
    else:
        print("✅ Network structure seems suitable for standard motif analysis")
    
    # Visualization suggestions
    print("\n=== VISUALIZATION SUGGESTIONS ===")
    if len(components) <= 20:
        print("Small number of components - can visualize full network")
    else:
        print("Many components - visualize largest few components separately")
    
    return {
        'components': components,
        'component_sizes': component_sizes,
        'densities': densities,
        'clustering': avg_clustering,
        'bias_indicators': bias_indicators
    }

def plot_component_analysis(component_sizes, densities):
    """Create plots to visualize component structure"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Component size distribution
    ax1.hist(component_sizes, bins=min(50, len(set(component_sizes))), alpha=0.7)
    ax1.set_xlabel('Component Size')
    ax1.set_ylabel('Number of Components')
    ax1.set_title('Distribution of Component Sizes')
    ax1.set_yscale('log')
    
    # Component size vs rank
    sorted_sizes = sorted(component_sizes, reverse=True)
    ax2.plot(range(1, len(sorted_sizes)+1), sorted_sizes, 'bo-', alpha=0.6)
    ax2.set_xlabel('Component Rank')
    ax2.set_ylabel('Component Size')
    ax2.set_title('Component Sizes (Ranked)')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    
    # Density distribution
    if densities:
        ax3.hist(densities, bins=20, alpha=0.7)
        ax3.set_xlabel('Density')
        ax3.set_ylabel('Number of Components')
        ax3.set_title('Distribution of Component Densities')
    
    # Size vs Density scatter
    if densities and len(densities) > 1:
        large_comp_sizes = [len(c) for c in component_sizes if len(c) >= 4][:len(densities)]
        ax4.scatter(large_comp_sizes, densities, alpha=0.6)
        ax4.set_xlabel('Component Size')
        ax4.set_ylabel('Component Density')
        ax4.set_title('Component Size vs Density')
    
    plt.tight_layout()
    plt.savefig('network_component_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python network_analysis.py <edgelist_file>")
        print("Example: python network_analysis.py loop_network.txt")
        sys.exit(1)
    
    edgelist_file = sys.argv[1]
    
    try:
        results = analyze_network_structure(edgelist_file)
        
        # Create plots
        plot_component_analysis(results['component_sizes'], results['densities'])
        
        print(f"\nAnalysis complete! Plots saved as 'network_component_analysis.png'")
        
    except Exception as e:
        print(f"Error analyzing network: {e}")
        print("Make sure your file is in the correct format (space/tab separated: source target)")