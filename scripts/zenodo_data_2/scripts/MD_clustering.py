import pandas as pd
import numpy as np
import copy
import os



def cluster_trajectory(replica_path):
    #
    # Calculate the cluster
    #
    from MDAnalysis.analysis import encore
    from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm
    import MDAnalysis as mda

    print()
    print()
    print("Start clustering of:", replica_path)


    u = mda.Universe(replica_path + "/cyp17a1.pdb", replica_path + "/cyp17a1.xtc")
    print("The sliced trajectory contains",len(u.trajectory),"frames")
    
    clustering_method = clm.AffinityPropagationNative(preference=-50.0,
                                                damping=0.9,
                                                max_iter=500,
                                                convergence_iter=50)
    print ('Clustering begins!')
    cluster_collection = encore.cluster([u], method=clustering_method)
    print(len(cluster_collection.clusters), "clusters have been found.")

    info_path = os.path.join(replica_path, "clusters_info.txt")
    with mda.Writer(replica_path + "clusters_cyp17a1.xtc", u.atoms.n_atoms) as w, \
         open(info_path, "w") as info_file:
        
        info_file.write("ClusterNumber\tCentroidFrameID\tClusterSize\n")  # Header

        for c_count, cluster in enumerate(cluster_collection.clusters):
            cluster_weight = cluster.size
            centroid_frameID = cluster.centroid
            centroid_frameNr = 1 + u.trajectory[cluster.centroid].frame
            print("Cluster found:", centroid_frameNr)
                   
            # Save the centroid
            u.trajectory[centroid_frameID]
            w.write(u.atoms)

            # Write info to file
            info_file.write(f"{c_count + 1}\t{centroid_frameNr}\t{cluster_weight}\n")
    
    print("Clustering done!")      
    print()

    # Print contents of the info file
    print("Contents of clusters_info.txt:")
    with open(info_path, "r") as f:
        print(f.read())

    print()
    print()  
    return 0
    

replica_path = "./"
cluster_trajectory(replica_path)

