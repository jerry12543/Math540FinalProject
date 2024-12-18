import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
import pandas as pd
import pathlib


def subset_expression_data(expression_data, sample_ids):
    return expression_data[[col for col in expression_data.columns if col in sample_ids]]

def build_network(expression_data, threshold=0.3):
    transposed_data = expression_data.T
    correlation_matrix = transposed_data.corr()
    adjacency_matrix = (correlation_matrix.abs() > threshold).astype(float) * correlation_matrix
    sparse_matrix = csr_matrix(adjacency_matrix.values)
    graph = nx.from_scipy_sparse_array(sparse_matrix)
    graph = nx.relabel_nodes(graph, {i: gene for i, gene in enumerate(expression_data.index)})
    return graph

def visualize_heatmap(graph, title="Heat Map of Network", save_path=None):
    adjacency_matrix = nx.to_numpy_array(graph)
    plt.figure(figsize=(10, 8))
    sns.heatmap(adjacency_matrix, cmap="viridis", cbar=True, square=True)
    plt.title(title)
    plt.xlabel("Genes")
    plt.ylabel("Genes")
    if save_path:
        pathlib.Path(save_path).parent.mkdir(parents=True, exist_ok=True) 
        plt.savefig(save_path, format="png", dpi=300)
        print(f"Heatmap saved to {save_path}")

def calculate_variability(expression_data, annotations, n_genes=500, gene_type_filter="protein-coding"):
    expression_data["Variance"] = expression_data.var(axis=1)
    merged_data = pd.merge(annotations, expression_data, left_index=True, right_index=True, how="inner")
    if gene_type_filter:
        merged_data = merged_data[merged_data["GeneType"] == gene_type_filter]
    most_variable_genes = merged_data.nlargest(n_genes, "Variance")
    return most_variable_genes[expression_data.columns].drop(columns=["Variance"])

def calc_network_stats(graph, annotations):
    degree_centrality = nx.degree_centrality(graph)
    betweenness_centrality = nx.betweenness_centrality(graph)
    closeness_centrality = nx.closeness_centrality(graph)
    eigenvector_centrality = nx.eigenvector_centrality(graph, max_iter=1000)

    top_nodes = sorted(degree_centrality, key=degree_centrality.get, reverse=True)[:10]
    results = []
    for node in top_nodes:
        symbol = annotations.loc[node, "Symbol"] if node in annotations.index else "N/A"
        results.append({
            "GeneID": node,
            "Symbol": symbol,
            "Degree": degree_centrality[node],
            "Betweenness": betweenness_centrality[node],
            "Closeness": closeness_centrality[node],
            "Eigenvector": eigenvector_centrality[node]
        })

    clustering_coefficients = nx.clustering(graph)
    avg_clustering = np.mean(list(clustering_coefficients.values()))
    return results, avg_clustering

def global_network_stats(graph):
    stats = {
        'Nodes': graph.number_of_nodes(),
        'Edges': graph.number_of_edges(),
        'Avg Degree': np.mean([deg for _, deg in graph.degree()]),
        'Density': nx.density(graph),
        'Avg Clustering': nx.average_clustering(graph),
        'Avg Shortest Path': nx.average_shortest_path_length(graph) if nx.is_connected(graph) else "N/A"
    }
    return stats

def compare_networks(grouped_samples, tpm_counts, annotations):
    network_metrics = {}

    for key, sample_ids in grouped_samples.items():
        filtered_data = subset_expression_data(tpm_counts, sample_ids)
        graph = build_network(filtered_data)
        stats = global_network_stats(graph)
        network_metrics[key] = stats

        centrality_results, avg_clustering = calc_network_stats(graph, annotations)
        print_latex_table(centrality_results, avg_clustering, f"{key} Network Centralities")

    # Print LaTeX Table
    print("\\begin{table}[h!]")
    print("\\centering")
    print("\\begin{tabular}{|l|c|c|c|c|c|c|}")
    print("\\hline")
    print("Condition & Nodes & Edges & Avg Degree & Density & Clustering & Shortest Path \\\\ \\hline")
    for key, stats in network_metrics.items():
        print(f"{key} & {stats['Nodes']} & {stats['Edges']} & {stats['Avg Degree']:.2f} & "
              f"{stats['Density']:.4f} & {stats['Avg Clustering']:.4f} & {stats['Avg Shortest Path']} \\\\ \\hline")
    print("\\end{tabular}")
    print("\\caption{Comparison of Network Topologies}")
    print("\\end{table}")

def print_latex_table(centrality_results, avg_clustering, title):
    print(f"\\section*{{{title}}}")
    print("\\begin{table}[h!]")
    print("\\centering")
    print("\\begin{tabular}{|l|l|c|c|c|c|}")
    print("\\hline")
    print("GeneID & Symbol & Degree & Betweenness & Closeness & Eigenvector \\\\ \\hline")
    for row in centrality_results:
        print(f"{row['GeneID']} & {row['Symbol']} & {row['Degree']:.4f} & {row['Betweenness']:.4f} & "
              f"{row['Closeness']:.4f} & {row['Eigenvector']:.4f} \\\\ \\hline")
    print("\\end{tabular}")
    print(f"\\caption{{Top Genes with Average Clustering Coefficient: {avg_clustering:.4f}}}")
    print("\\end{table}\n")


def analyze_degree_centrality(graph, annotations, top_n=10):
    """
    Analyze degree centrality and identify hub genes.
    
    Args:
        graph (nx.Graph): Network graph.
        annotations (pd.DataFrame): Gene annotations.
        top_n (int): Number of top hub genes to return.
    
    Returns:
        pd.DataFrame: Dataframe of top hub genes with centrality values.
    """
    # Calculate degree centrality
    degree_centrality = nx.degree_centrality(graph)

    # Sort and select top N genes
    sorted_genes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    # Map gene symbols and centrality scores
    hub_genes = pd.DataFrame([
        {"GeneID": gene, 
         "Symbol": annotations.loc[gene, "Symbol"] if gene in annotations.index else "N/A", 
         "DegreeCentrality": centrality}
        for gene, centrality in sorted_genes
    ])
    return hub_genes

def plot_degree_centrality_distribution(graph, title, save_path=None):
    """
    Plot the degree centrality distribution as a histogram.
    
    Args:
        graph (nx.Graph): Network graph.
        title (str): Title for the plot.
        save_path (str): Path to save the plot (optional).
    """
    degree_centrality = nx.degree_centrality(graph).values()
    plt.figure(figsize=(8, 6))
    plt.hist(degree_centrality, bins=50, color='blue', alpha=0.7)
    plt.title(f"Degree Centrality Distribution: {title}")
    plt.xlabel("Degree Centrality")
    plt.ylabel("Frequency")
    if save_path:
        pathlib.Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300)
        print(f"Saved histogram to {save_path}")
    plt.show()
