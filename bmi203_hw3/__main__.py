import sys
from .io import read_active_sites, write_H_clustering, write_P_clustering
from .cluster import cluster_by_partitioning, cluster_hierarchically, silhouette, create_similarity_matrix

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    sim_matrix = create_similarity_matrix(active_sites)
    clustering = cluster_by_partitioning(active_sites,sim_matrix)
    print clustering
    print type(clustering)
    print type(clustering[0])
    s = silhouette(clustering,active_sites,sim_matrix)
    write_P_clustering(sys.argv[3], clustering, s)

if sys.argv[1][0:2] == '-H':
    print("Clustering using Hierarchical method")
    sim_matrix = create_similarity_matrix(active_sites)
    clustering = cluster_hierarchically(active_sites,sim_matrix)
    write_H_clustering(sys.argv[3], clustering)
