'''
    similarity_matrix = -pairwise_distances(feature_vectors, metric='euclidean').astype(np.float64)

    # Compare with your precomputed matrix
    print("Difference between computed and precomputed:", np.sum(computed_similarity - precomputed_similarity_matrix))

    # Perform Affinity Propagation using a precomputed similarity matrix
    clustering = AffinityPropagation(affinity='precomputed', preference=-4000, damping=0.5, max_iter=200, convergence_iter=15)
    clustering.fit(similarity_matrix)
'''


import sys
import numpy as np
import os
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import pairwise_distances

def perform_clustering(input_file, output_file):
    all_data = []
    file_sources = []
    contributing_snapshots = {}
    
    # Read data from the input file
    with open(input_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) > 3:
                snapshot_center = int(parts[0])
                snapshot_name = parts[1]
                num_snapshots = int(parts[2])  # Number of contributing snapshots
                feature_vector = list(map(float, parts[3:]))
                all_data.append((snapshot_center, snapshot_name, num_snapshots, feature_vector))
                file_sources.append(snapshot_name)
                if snapshot_name not in contributing_snapshots:
                    contributing_snapshots[snapshot_name] = {}
                contributing_snapshots[snapshot_name][snapshot_center] = num_snapshots
    
    structure_indices = [data[0] for data in all_data]
    feature_vectors = np.array([data[3] for data in all_data])
    file_names = [data[1] for data in all_data]

    
    print(f"Feature vector shape: {feature_vectors.shape}")

    #prije max_iter 200 con 15 
    #Perform Affinity Propagation clustering
    clustering = AffinityPropagation(preference=-22000, damping=0.5, max_iter=500, convergence_iter=100)
    clustering.fit(feature_vectors)
    
    cluster_centers_indices = clustering.cluster_centers_indices_
    labels = clustering.labels_
    
    cluster_counts = {}
    cluster_snapshots = {}
    
    for cluster_id in range(len(cluster_centers_indices)):
        cluster_counts[cluster_id] = {}
        cluster_snapshots[cluster_id] = {}
    
    for index, label in enumerate(labels):
        file_name = file_names[index]
        snapshot_id = structure_indices[index]
        num_snapshots = contributing_snapshots[file_name][snapshot_id]
        
        if file_name not in cluster_counts[label]:
            cluster_counts[label][file_name] = 0
        cluster_counts[label][file_name] += 1
        
        if file_name not in cluster_snapshots[label]:
            cluster_snapshots[label][file_name] = 0
        cluster_snapshots[label][file_name] += num_snapshots  # Track per file_name
    
    # Save clustering results to a file
    with open(output_file, "w") as f:
        f.write(f"Total clusters found: {len(cluster_centers_indices)}\n\n")

        for cluster_id, centroid_index in enumerate(cluster_centers_indices):
            cluster_members = np.where(labels == cluster_id)[0]
            cluster_size = len(cluster_members)
            
            f.write(f"Cluster {cluster_id}:\n")
            f.write(f"  Centroid Structure Index: {structure_indices[centroid_index]}\n")
            f.write(f"  Centroid File Origin: {file_names[centroid_index]}\n")
            f.write(f"  Number of Structures: {cluster_size}\n")
            
            f.write("  File Distribution: \n")
            for file_name, count in cluster_counts[cluster_id].items():
                f.write(f"    {file_name}: {count} (Contributing Snapshots: {cluster_snapshots[cluster_id][file_name]})\n")
            
            f.write(f"  Cluster Members: {', '.join(map(str, [structure_indices[i] for i in cluster_members]))}\n")
            f.write("\n")
    
    print(f"Clustering complete. Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clustering_script.py cluster_vectors.txt output.txt")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        perform_clustering(input_file, output_file)

