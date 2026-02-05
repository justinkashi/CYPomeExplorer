import os
import numpy as np

def extract_centroids(cluster_filename):
    """Extract centroid indices from a cluster file using a generator."""
    with open(cluster_filename, 'r') as file:
        for line in file:
            if line.startswith("  Centroid Structure Index: "):
                yield line.split(": ")[1].strip()
def extract_snapshot_number(cluster_filename):
    '''Extracts the number of snapshots populating each cluster'''
    with open(cluster_filename,'r') as file:
        numbers=[]
        for line in file:
            if line.startswith("  Number of Structures:"):
                numbers.append(line.split(": ")[1].strip())
    return (np.array(numbers))
def read_vectors(cluster_filename):
    """Read all vector lengths for centroid structures from a file and return as a list of NumPy arrays."""
    vectors_list = []
    with open(cluster_filename, 'r') as file:
        for line in file:
            if line.strip().startswith("Vectors lengths for centroid structure:"):
                values_str = line.split(":", 1)[1].strip()

                while not values_str.endswith("]"):
                    next_line = next(file).strip()
                    values_str += " " + next_line

                values = eval(values_str)  # controlled input, safe use here
                vectors_list.append(np.array(values))

    return vectors_list


def process_cluster_files(clusters, output_filename):
    """Process cluster files and write matching centroids to output."""
    results = []
    for cluster_file in clusters:
        vectors = read_vectors(cluster_file)
        cluster_filename = os.path.basename(cluster_file)  # Get filename only
        numbers=extract_snapshot_number(cluster_file)
        i=0
        centroids = extract_centroids(cluster_file)
        for centroid in centroids:
            vector_str = ' '.join(f"{x:.3f}" for x in vectors[i])  # format floats nicely
            results.append(f"{centroid} {cluster_filename} {numbers[i]} {vector_str}\n")
            i+=1

    # Write everything at once to improve performance
    with open(output_filename, 'w') as output_file:
        output_file.writelines(results)

if __name__ == "__main__":
    clusters = [
 	        ("clusters_3a4.txt"),
       ("clusters_1a2.txt"),
        ("clusters_2d6.txt"),
        ("clusters_2c9.txt"),
        ("clusters_2c19.txt"),
        
        ("clusters_81a1.txt"),
        ("clusters_81a2.txt"),
        ("clusters_81a4.txt"),
        ("clusters_81a9.txt"),
        ("clusters_81a16.txt"),
        
        ("clusters_72a31.txt"),
        ("clusters_81a6.txt"),
        
        ("clusters_72a188.txt"),
        ("clusters_72a208.txt"),
        ("clusters_79a1.txt"),
        ("clusters_79e1.txt"),
        ("clusters_81f2.txt"),
        ("clusters_81f4.txt"),
        ("clusters_90c1.txt"),
        ("clusters_90d1.txt"),
        ("clusters_2a6.txt")
    ]  # Add all corresponding cluster and vector file pairs here
    output_file = "clusters_all_vectors.txt"
    process_cluster_files(clusters, output_file)

