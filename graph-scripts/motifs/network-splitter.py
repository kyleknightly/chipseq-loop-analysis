import pickle
import networkx as nx
from collections import defaultdict
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
import community as community_louvain  # python-louvain package
from typing import Dict, List, Tuple, Set
import os
from pathlib import Path

class ChromosomeGraphSplitter:
    """
    A class to split genomic position graphs by chromosome and then optimally
    partition each chromosome into components while minimizing inter-component edges.
    """
    
    def __init__(self, graph_path: str):
        """Load the pickled graph and create integer node mapping."""
        with open(graph_path, 'rb') as f:
            original_graph = pickle.load(f)
        
        # Validate that nodes are tuples of (chromosome, start, end)
        sample_node = next(iter(original_graph.nodes()))
        if not isinstance(sample_node, tuple) or len(sample_node) != 3:
            raise ValueError("Expected nodes to be tuples of (chromosome, start, end)")
        
        print("Converting tuple node names to integer IDs...")
        
        # Create mappings between original tuple nodes and integer IDs
        self.original_nodes = list(original_graph.nodes())
        self.node_to_id = {node: i for i, node in enumerate(self.original_nodes)}
        self.id_to_node = {i: node for i, node in enumerate(self.original_nodes)}
        
        # Create new graph with integer node IDs
        self.graph = nx.Graph()
        
        # Add nodes with integer IDs
        self.graph.add_nodes_from(range(len(self.original_nodes)))
        
        # Add edges with converted IDs
        for edge in original_graph.edges():
            id1 = self.node_to_id[edge[0]]
            id2 = self.node_to_id[edge[1]]
            self.graph.add_edge(id1, id2)
        
        # Create chromosome mapping for efficient splitting
        self.node_chromosomes = {}
        self.chromosome_nodes = defaultdict(list)
        
        for node_id, original_node in self.id_to_node.items():
            chromosome = original_node[0]
            self.node_chromosomes[node_id] = chromosome
            self.chromosome_nodes[chromosome].append(node_id)
        
        print(f"Converted {len(self.original_nodes)} tuple nodes to integer IDs")
        print(f"Found {len(self.chromosome_nodes)} chromosomes")
    
    def get_original_node(self, node_id: int) -> Tuple[str, int, int]:
        """Convert integer node ID back to original tuple."""
        return self.id_to_node[node_id]
    
    def get_node_id(self, original_node: Tuple[str, int, int]) -> int:
        """Convert original tuple node to integer ID."""
        return self.node_to_id[original_node]
    
    def split_by_chromosome(self) -> Dict[str, nx.Graph]:
        """
        Split the main graph into separate subgraphs for each chromosome.
        
        Returns:
            Dictionary mapping chromosome names to their respective subgraphs
        """
        subgraphs = {}
        
        for chromosome, node_ids in self.chromosome_nodes.items():
            subgraph = self.graph.subgraph(node_ids).copy()
            subgraphs[chromosome] = subgraph
            print(f"Chromosome {chromosome}: {len(node_ids)} nodes, {subgraph.number_of_edges()} edges")
        
        return subgraphs
    
    def evaluate_partition_quality(self, graph: nx.Graph, partition: Dict) -> Tuple[float, int]:
        """
        Evaluate the quality of a partition using modularity and edge cut metrics.
        
        Args:
            graph: The graph being partitioned
            partition: Dictionary mapping nodes to community/cluster IDs
            
        Returns:
            Tuple of (modularity, edges_between_communities)
        """
        # Calculate modularity
        modularity = community_louvain.modularity(partition, graph)
        
        # Count edges between different communities
        edges_between = 0
        for edge in graph.edges():
            if partition[edge[0]] != partition[edge[1]]:
                edges_between += 1
        
        return modularity, edges_between
    
    def louvain_clustering(self, graph: nx.Graph, resolution: float = 1.0) -> Dict:
        """
        Apply Louvain community detection algorithm.
        
        Args:
            graph: Input graph
            resolution: Resolution parameter for community detection
            
        Returns:
            Dictionary mapping nodes to community IDs
        """
        if graph.number_of_nodes() == 0:
            return {}
        
        partition = community_louvain.best_partition(graph, resolution=resolution, random_state=42)
        return partition
    
    def spectral_clustering_with_optimization(self, graph: nx.Graph, max_clusters: int = None) -> Dict:
        """
        Apply spectral clustering with automatic cluster number selection.
        
        Args:
            graph: Input graph
            max_clusters: Maximum number of clusters to consider
            
        Returns:
            Dictionary mapping nodes to cluster IDs
        """
        if graph.number_of_nodes() == 0:
            return {}
        
        if graph.number_of_nodes() == 1:
            return {list(graph.nodes())[0]: 0}
        
        # Convert to adjacency matrix
        adj_matrix = nx.adjacency_matrix(graph).toarray()
        
        if max_clusters is None:
            max_clusters = min(20, graph.number_of_nodes() // 2)
        
        best_score = -1
        best_partition = None
        best_n_clusters = 2
        
        # Try different numbers of clusters
        for n_clusters in range(2, min(max_clusters + 1, graph.number_of_nodes())):
            try:
                clustering = SpectralClustering(
                    n_clusters=n_clusters, 
                    affinity='precomputed',
                    random_state=42,
                    assign_labels='kmeans'
                )
                
                cluster_labels = clustering.fit_predict(adj_matrix)
                
                # Calculate silhouette score
                if len(set(cluster_labels)) > 1:
                    score = silhouette_score(adj_matrix, cluster_labels, metric='precomputed')
                    
                    if score > best_score:
                        best_score = score
                        best_n_clusters = n_clusters
                        best_partition = {node: int(label) for node, label in 
                                        zip(graph.nodes(), cluster_labels)}
                        
            except Exception as e:
                print(f"Error with {n_clusters} clusters: {e}")
                continue
        
        if best_partition is None:
            # Fallback: assign all nodes to one cluster
            best_partition = {node: 0 for node in graph.nodes()}
        
        return best_partition
    
    def optimize_resolution_louvain(self, graph: nx.Graph, 
                                  resolution_range: Tuple[float, float] = (0.5, 2.0),
                                  num_trials: int = 20) -> Dict:
        """
        Find optimal resolution parameter for Louvain algorithm.
        
        Args:
            graph: Input graph
            resolution_range: Range of resolution values to try
            num_trials: Number of resolution values to test
            
        Returns:
            Best partition found
        """
        if graph.number_of_nodes() == 0:
            return {}
        
        resolutions = np.linspace(resolution_range[0], resolution_range[1], num_trials)
        best_modularity = -1
        best_partition = None
        
        for resolution in resolutions:
            partition = self.louvain_clustering(graph, resolution)
            modularity, edges_between = self.evaluate_partition_quality(graph, partition)
            
            if modularity > best_modularity:
                best_modularity = modularity
                best_partition = partition
        
        return best_partition
    
    def partition_chromosome_graph(self, graph: nx.Graph, method: str = 'louvain_optimized') -> Dict:
        """
        Partition a chromosome graph into optimal components.
        
        Args:
            graph: Chromosome subgraph
            method: Partitioning method ('louvain', 'louvain_optimized', 'spectral')
            
        Returns:
            Dictionary mapping nodes to component IDs
        """
        if graph.number_of_nodes() == 0:
            return {}
        
        if method == 'louvain':
            return self.louvain_clustering(graph)
        elif method == 'louvain_optimized':
            return self.optimize_resolution_louvain(graph)
        elif method == 'spectral':
            return self.spectral_clustering_with_optimization(graph)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def split_and_partition_all(self, method: str = 'louvain_optimized') -> Dict[str, Dict]:
        """
        Split by chromosome and partition each chromosome optimally.
        
        Args:
            method: Partitioning method to use
            
        Returns:
            Dictionary mapping chromosome -> partition dictionary
        """
        print("Splitting graph by chromosome...")
        chromosome_graphs = self.split_by_chromosome()
        
        results = {}
        print(f"\nPartitioning each chromosome using {method} method...")
        
        for chromosome, subgraph in chromosome_graphs.items():
            print(f"\nProcessing chromosome {chromosome}...")
            
            if subgraph.number_of_nodes() == 0:
                results[chromosome] = {}
                continue
            
            partition = self.partition_chromosome_graph(subgraph, method)
            modularity, edges_between = self.evaluate_partition_quality(subgraph, partition)
            
            num_components = len(set(partition.values()))
            print(f"  - Found {num_components} components")
            print(f"  - Modularity: {modularity:.4f}")
            print(f"  - Edges between components: {edges_between}")
            print(f"  - Edge cut ratio: {edges_between/subgraph.number_of_edges():.4f}")
            
            results[chromosome] = partition
        
        return results
    
    def extract_components_as_graphs(self, chromosome_partitions: Dict[str, Dict]) -> Dict[str, List[nx.Graph]]:
        """
        Extract each component as a separate graph.
        
        Args:
            chromosome_partitions: Result from split_and_partition_all()
            
        Returns:
            Dictionary mapping chromosome -> list of component graphs
        """
        chromosome_graphs = self.split_by_chromosome()
        component_graphs = {}
        
        for chromosome, partition in chromosome_partitions.items():
            if not partition:
                component_graphs[chromosome] = []
                continue
            
            subgraph = chromosome_graphs[chromosome]
            
            # Group nodes by component
            components = defaultdict(list)
            for node, component_id in partition.items():
                components[component_id].append(node)
            
            # Create subgraph for each component
            graphs = []
            for component_id, nodes in components.items():
                component_graph = subgraph.subgraph(nodes).copy()
                graphs.append(component_graph)
            
            component_graphs[chromosome] = graphs
        
        return component_graphs
    
    def save_results(self, chromosome_partitions: Dict[str, Dict], 
                    outdir: str, filename_prefix: str = 'chromosome_partition'):
        """
        Save partitioning results to specified output directory in mfinder format.
        
        Args:
            chromosome_partitions: Result from split_and_partition_all()
            outdir: Output directory path
            filename_prefix: Prefix for output files
        """
        # Create output directory if it doesn't exist
        Path(outdir).mkdir(parents=True, exist_ok=True)
        
        # Save partition mappings as pickle (for reference)
        mappings_path = os.path.join(outdir, f'{filename_prefix}_mappings.pkl')
        with open(mappings_path, 'wb') as f:
            pickle.dump(chromosome_partitions, f)
        
        # Save chromosome graphs in mfinder format
        chromosome_graphs = self.split_by_chromosome()
        chr_dir = os.path.join(outdir, 'chromosome_networks')
        Path(chr_dir).mkdir(parents=True, exist_ok=True)
        
        for chromosome, graph in chromosome_graphs.items():
            chr_path = os.path.join(chr_dir, f'{chromosome}.txt')
            self._save_graph_mfinder_format(graph, chr_path)
        
        # Save component graphs in mfinder format
        component_graphs = self.extract_components_as_graphs(chromosome_partitions)
        components_dir = os.path.join(outdir, 'component_networks')
        Path(components_dir).mkdir(parents=True, exist_ok=True)
        
        for chromosome, graphs in component_graphs.items():
            for i, component_graph in enumerate(graphs):
                comp_path = os.path.join(components_dir, f'{chromosome}_component_{i}.txt')
                self._save_graph_mfinder_format(component_graph, comp_path)
        
        # Save node mappings (original nodes to integer IDs used in mfinder files)
        self._save_node_mappings(chromosome_partitions, outdir, filename_prefix)
        
        # Save summary statistics
        self._save_summary_stats(chromosome_partitions, outdir, filename_prefix)
        
        print(f"Results saved to {outdir}:")
        print(f"  - Partition mappings: {mappings_path}")
        print(f"  - Chromosome networks (mfinder format): {chr_dir}/")
        print(f"  - Component networks (mfinder format): {components_dir}/")
        print(f"  - Node mappings: {os.path.join(outdir, f'{filename_prefix}_node_mappings.txt')}")
        print(f"  - Summary statistics: {os.path.join(outdir, f'{filename_prefix}_summary.txt')}")
    
    def _save_graph_mfinder_format(self, graph: nx.Graph, filepath: str):
        """
        Save a graph in mfinder format (tab-separated edges).
        
        Args:
            graph: NetworkX graph with integer node IDs to save
            filepath: Output file path
        """
        if graph.number_of_nodes() == 0:
            # Create empty file for empty graphs
            with open(filepath, 'w') as f:
                f.write("# Empty graph\n")
            return
        
        # Renumber nodes to be consecutive starting from 1 for mfinder
        nodes = list(graph.nodes())
        node_renumber = {node: i+1 for i, node in enumerate(nodes)}
        
        with open(filepath, 'w') as f:
            # Write header comment
            f.write(f"# Graph with {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges\n")
            f.write(f"# Edges in format: node1_id<tab>node2_id\n")
            
            # Write edges in tab-separated format
            for edge in graph.edges():
                node1_id = node_renumber[edge[0]]
                node2_id = node_renumber[edge[1]]
                f.write(f"{node1_id}\t{node2_id}\n")
    
    def _save_node_mappings(self, chromosome_partitions: Dict[str, Dict], 
                           outdir: str, filename_prefix: str):
        """
        Save comprehensive node mappings from mfinder IDs back to genomic positions.
        
        Args:
            chromosome_partitions: Result from split_and_partition_all()
            outdir: Output directory path
            filename_prefix: Prefix for output files
        """
        # Save master lookup table (internal ID -> genomic position)
        master_mapping_path = os.path.join(outdir, f'{filename_prefix}_master_node_lookup.txt')
        with open(master_mapping_path, 'w') as f:
            f.write("# Master node lookup table\n")
            f.write("# Format: internal_id<tab>chromosome<tab>start<tab>end\n")
            
            for internal_id, original_node in self.id_to_node.items():
                chr_name, start, end = original_node
                f.write(f"{internal_id}\t{chr_name}\t{start}\t{end}\n")
        
        chromosome_graphs = self.split_by_chromosome()
        component_graphs = self.extract_components_as_graphs(chromosome_partitions)
        
        # Save chromosome-level mappings (mfinder ID -> internal ID -> genomic position)
        chr_mapping_path = os.path.join(outdir, f'{filename_prefix}_chromosome_mfinder_mappings.txt')
        with open(chr_mapping_path, 'w') as f:
            f.write("# Chromosome-level mfinder node mappings\n")
            f.write("# Format: chromosome<tab>mfinder_id<tab>internal_id<tab>chromosome<tab>start<tab>end\n")
            
            for chromosome, graph in chromosome_graphs.items():
                nodes = list(graph.nodes())
                mfinder_renumber = {node: i+1 for i, node in enumerate(nodes)}
                
                for internal_id, mfinder_id in mfinder_renumber.items():
                    original_node = self.get_original_node(internal_id)
                    chr_name, start, end = original_node
                    f.write(f"{chromosome}\t{mfinder_id}\t{internal_id}\t{chr_name}\t{start}\t{end}\n")
        
        # Save component-level mappings
        comp_mapping_path = os.path.join(outdir, f'{filename_prefix}_component_mfinder_mappings.txt')
        with open(comp_mapping_path, 'w') as f:
            f.write("# Component-level mfinder node mappings\n")
            f.write("# Format: chromosome<tab>component_id<tab>mfinder_id<tab>internal_id<tab>chromosome<tab>start<tab>end\n")
            
            for chromosome, graphs in component_graphs.items():
                for comp_id, graph in enumerate(graphs):
                    nodes = list(graph.nodes())
                    mfinder_renumber = {node: i+1 for i, node in enumerate(nodes)}
                    
                    for internal_id, mfinder_id in mfinder_renumber.items():
                        original_node = self.get_original_node(internal_id)
                        chr_name, start, end = original_node
                        f.write(f"{chromosome}\t{comp_id}\t{mfinder_id}\t{internal_id}\t{chr_name}\t{start}\t{end}\n")
    
    def _save_summary_stats(self, chromosome_partitions: Dict[str, Dict], 
                           outdir: str, filename_prefix: str):
        """Save summary statistics to a text file."""
        chromosome_graphs = self.split_by_chromosome()
        
        summary_path = os.path.join(outdir, f'{filename_prefix}_summary.txt')
        with open(summary_path, 'w') as f:
            f.write("Chromosome Graph Partitioning Summary\n")
            f.write("=" * 40 + "\n\n")
            
            total_nodes = sum(g.number_of_nodes() for g in chromosome_graphs.values())
            total_edges = sum(g.number_of_edges() for g in chromosome_graphs.values())
            total_components = sum(len(set(p.values())) for p in chromosome_partitions.values() if p)
            
            f.write(f"Overall Statistics:\n")
            f.write(f"  Total nodes: {total_nodes}\n")
            f.write(f"  Total edges: {total_edges}\n")
            f.write(f"  Total chromosomes: {len(chromosome_graphs)}\n")
            f.write(f"  Total components: {total_components}\n\n")
            
            f.write("Per-Chromosome Statistics:\n")
            f.write("-" * 30 + "\n")
            
            for chromosome in sorted(chromosome_graphs.keys()):
                subgraph = chromosome_graphs[chromosome]
                partition = chromosome_partitions.get(chromosome, {})
                
                if not partition:
                    f.write(f"{chromosome}: Empty partition\n")
                    continue
                
                num_components = len(set(partition.values()))
                modularity, edges_between = self.evaluate_partition_quality(subgraph, partition)
                edge_cut_ratio = edges_between / subgraph.number_of_edges() if subgraph.number_of_edges() > 0 else 0
                
                f.write(f"{chromosome}:\n")
                f.write(f"  Nodes: {subgraph.number_of_nodes()}\n")
                f.write(f"  Edges: {subgraph.number_of_edges()}\n")
                f.write(f"  Components: {num_components}\n")
                f.write(f"  Modularity: {modularity:.4f}\n")
                f.write(f"  Inter-component edges: {edges_between}\n")
                f.write(f"  Edge cut ratio: {edge_cut_ratio:.4f}\n\n")


# Main execution function for any graph
def process_graph(graph_path: str, outdir: str, method: str = 'louvain_optimized'):
    """
    Process any chromosome position graph with optimal partitioning.
    
    Args:
        graph_path: Path to the input graph (.gpickle file)
        outdir: Output directory for all results
        method: Partitioning method ('louvain_optimized', 'spectral', or 'louvain')
    
    Returns:
        ChromosomeGraphSplitter instance and results
    """
    print(f"Loading graph from: {graph_path}")
    splitter = ChromosomeGraphSplitter(graph_path)
    
    print(f"Total nodes: {splitter.graph.number_of_nodes()}")
    print(f"Total edges: {splitter.graph.number_of_edges()}")
    
    # Extract base filename for output prefix
    base_filename = os.path.splitext(os.path.basename(graph_path))[0]
    output_prefix = f'{base_filename}_{method}'
    
    # Perform the splitting and partitioning
    print(f"\nProcessing with method: {method}")
    results = splitter.split_and_partition_all(method=method)
    
    # Save all results to the specified output directory
    splitter.save_results(results, outdir, output_prefix)
    
    return splitter, results


# Example usage
def example_usage():
    """Example of how to process any graph."""
    
    # Example with your HepG2 graph
    graph_path = '/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/hepg2-anchor-graph.gpickle'
    output_directory = '/path/to/your/output/directory'
    
    # Process with Louvain optimization (recommended)
    splitter, results = process_graph(graph_path, output_directory, method='louvain_optimized')
    
    # Optionally compare with spectral clustering
    print("\n" + "="*50)
    print("Comparing with spectral clustering...")
    splitter_spectral, results_spectral = process_graph(
        graph_path, output_directory, method='spectral'
    )
    
    return splitter, results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Split and partition chromosome graphs')
    parser.add_argument('graph_path', help='Path to input graph (.gpickle file)')
    parser.add_argument('outdir', help='Output directory for results')
    parser.add_argument('--method', default='louvain_optimized', 
                       choices=['louvain', 'louvain_optimized', 'spectral'],
                       help='Partitioning method (default: louvain_optimized)')
    
    args = parser.parse_args()
    
    print(f"Processing graph: {args.graph_path}")
    print(f"Output directory: {args.outdir}")
    print(f"Method: {args.method}")
    
    splitter, results = process_graph(args.graph_path, args.outdir, args.method)
    
    print("\nProcessing complete! Check the output directory for results.")