#!/usr/bin/env python3
"""
Filter graph to extract largest connected components for motif analysis.
Usage: python filter_components.py input.gpickle [min_size]
"""

import networkx as nx
import sys
import os

def filter_components(input_file, min_size=100, save_largest_only=True):
    """
    Filter graph to extract meaningful components for motif analysis.
    
    Args:
        input_file: Path to .gpickle file
        min_size: Minimum component size to keep
        save_largest_only: If True, save only the largest component
    """
    
    print(f"Loading graph from {input_file}...")
    G = nx.read_gpickle(input_file)
    
    print(f"Original graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"Directed: {G.is_directed()}")
    
    # Get connected components (use appropriate function for directed/undirected)
    if G.is_directed():
        print("Analyzing weakly connected components...")
        components = list(nx.weakly_connected_components(G))
    else:
        print("Analyzing connected components...")
        components = list(nx.connected_components(G))
    
    # Sort components by size (largest first)
    components = sorted(components, key=len, reverse=True)
    component_sizes = [len(c) for c in components]
    
    print(f"\nComponent analysis:")
    print(f"Total components: {len(components)}")
    print(f"Largest component: {component_sizes[0]} nodes ({component_sizes[0]/G.number_of_nodes()*100:.1f}%)")
    print(f"Components with ≥{min_size} nodes: {sum(1 for s in component_sizes if s >= min_size)}")
    
    # Show top 10 component sizes
    print(f"\nTop 10 component sizes:")
    for i, size in enumerate(component_sizes[:10]):
        print(f"  Component {i+1}: {size} nodes")
    
    # Filter components
    if save_largest_only:
        # Save only the largest component
        largest_component = components[0]
        largest_subgraph = G.subgraph(largest_component)
        
        output_base = os.path.splitext(input_file)[0]
        output_file = f"{output_base}_LCC.txt"
        
        # Save as edgelist for mfinder
        nx.write_edgelist(largest_subgraph, output_file, data=False)
        
        print(f"\n✅ Saved largest component:")
        print(f"   File: {output_file}")
        print(f"   Size: {largest_subgraph.number_of_nodes()} nodes, {largest_subgraph.number_of_edges()} edges")
        print(f"   Density: {nx.density(largest_subgraph):.6f}")
        
        # Also save as gpickle for further analysis (create new graph to avoid pickle issues)
        try:
            gpickle_output = f"{output_base}_LCC.gpickle"
            # Create a new graph instead of using subgraph to avoid pickle issues
            new_graph = nx.DiGraph() if G.is_directed() else nx.Graph()
            new_graph.add_nodes_from(largest_subgraph.nodes(data=True))
            new_graph.add_edges_from(largest_subgraph.edges(data=True))
            nx.write_gpickle(new_graph, gpickle_output)
            print(f"   Also saved as: {gpickle_output}")
        except Exception as e:
            print(f"   (Couldn't save gpickle: {e})")
        
        return [(output_file, largest_subgraph)]
    
    else:
        # Save all components above min_size
        large_components = [c for c in components if len(c) >= min_size]
        
        output_base = os.path.splitext(input_file)[0]
        saved_files = []
        
        print(f"\n✅ Saving {len(large_components)} components with ≥{min_size} nodes:")
        
        for i, component in enumerate(large_components):
            subgraph = G.subgraph(component)
            output_file = f"{output_base}_component_{i+1}.txt"
            
            # Save as edgelist for mfinder
            nx.write_edgelist(subgraph, output_file, data=False)
            
            print(f"   Component {i+1}: {output_file}")
            print(f"     Size: {subgraph.number_of_nodes()} nodes, {subgraph.number_of_edges()} edges")
            print(f"     Density: {nx.density(subgraph):.6f}")
            
            saved_files.append((output_file, subgraph))
        
        return saved_files

def main():
    if len(sys.argv) < 2:
        print("Usage: python filter_components.py input.gpickle [min_size] [--all]")
        print("  input.gpickle: Input graph file")
        print("  min_size: Minimum component size (default: 100)")
        print("  --all: Save all large components instead of just largest")
        print("\nExample:")
        print("  python filter_components.py graph.gpickle")
        print("  python filter_components.py graph.gpickle 50 --all")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Parse arguments
    min_size = 100
    save_largest_only = True
    
    for arg in sys.argv[2:]:
        if arg == "--all":
            save_largest_only = False
        elif arg.isdigit():
            min_size = int(arg)
    
    # Validate input file
    if not os.path.exists(input_file):
        print(f"Error: File {input_file} not found")
        sys.exit(1)
    
    if not input_file.endswith(('.gpickle', '.pickle')):
        print(f"Warning: File doesn't end with .gpickle or .pickle")
    
    try:
        saved_files = filter_components(input_file, min_size, save_largest_only)
        
        print(f"\n🎯 Next steps for motif analysis:")
        for output_file, subgraph in saved_files:
            print(f"\n   mfinder {output_file} -s 4 -r 100 -nd")
        
        print(f"\n💡 Tips:")
        print(f"   - Start with the largest component first")
        print(f"   - Use more random networks (-r 1000) for better statistics")
        print(f"   - Compare results between components")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()