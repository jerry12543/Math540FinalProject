import loader
from collections import defaultdict
from tables import calculate_variability, compare_networks, analyze_degree_centrality, plot_degree_centrality_distribution, subset_expression_data, build_network
import config

# Load data
metadata = loader.metadata
raw_counts = loader.raw_counts
tpm_counts = loader.tpm_counts
annotations = loader.annotations

# Main Execution
most_variable_genes = calculate_variability(tpm_counts, annotations)
tpm_counts = most_variable_genes

grouped_samples = defaultdict(list)
for gsm_id, attributes in metadata.items():
    diagnosis = attributes.get("diagnosis")
    cell_type = attributes.get("cell_type")
    grouped_samples[(diagnosis, cell_type)].append(gsm_id)

# Process each group and analyze hubs
for key, sample_ids in grouped_samples.items():
    print(f"Processing {key}...")
    
    # Subset data and build network
    filtered_data = subset_expression_data(tpm_counts, sample_ids)
    graph = build_network(filtered_data)
    
    # Analyze hubs
    hub_genes = analyze_degree_centrality(graph, annotations)
    print(f"\nTop Hub Genes in {key}:")
    print(hub_genes.to_string(index=False))
    
    # Visualize degree centrality distribution
    save_path = config.HEATMAPSPATH / f"{key}_degree_centrality_hist.png"
    plot_degree_centrality_distribution(graph, title=f"{key} Network", save_path=save_path)

# Compare networks and generate results
compare_networks(grouped_samples, tpm_counts, annotations)

