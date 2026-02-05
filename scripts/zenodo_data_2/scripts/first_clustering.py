import sys
import numpy as np
from sklearn.cluster import AffinityPropagation

def perform_clustering(input_file, output_file):
    # Load the dataset
    with open(input_file, "r") as f:
        data = [list(map(float, line.strip().split())) for line in f]

    # Extract structure indices and feature vectors
    structure_indices = [int(row[0]) for row in data]
    feature_vectors = np.array([row[1:] for row in data])
    print(f"Feature vector shape: {feature_vectors.shape}")

    # Perform Affinity Propagation clustering
    clustering = AffinityPropagation(preference=-2000, damping=0.9, max_iter=500, convergence_iter=100)
    clustering.fit(feature_vectors)

    # Extract clustering results
    cluster_centers_indices = clustering.cluster_centers_indices_
    labels = clustering.labels_

    # Save clustering results to a file
    with open(output_file, "w") as f:
        f.write(f"Total clusters found: {len(cluster_centers_indices)}\n\n")

        for cluster_id, centroid_index in enumerate(cluster_centers_indices):
            # Find all frames assigned to this cluster
            cluster_members = np.where(labels == cluster_id)[0]
            cluster_size = len(cluster_members)  # Number of structures in the cluster

            # Extract feature vector of centroid structure
            centroid_features = feature_vectors[centroid_index]

            # Write cluster information
            f.write(f"Cluster {cluster_id}:\n")
            f.write(f"  Centroid Structure Index: {structure_indices[centroid_index]}\n")
            f.write(f"  Number of Structures: {cluster_size}\n")
            f.write(f"  Vectors lengths for centroid structure: {centroid_features.tolist()}\n\n")

    print(f"Clustering complete. Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clustering_script.py <input_file.txt> <output_file.txt>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        perform_clustering(input_file, output_file)


