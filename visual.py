from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
import numpy as np
import matplotlib.pyplot as plt

def visualize():
    """
    This function visualizes a dendogram showing the evolutinary tree adn the relationship between the number of clusters and Silhouette score
    """

    # Cluster using mean linkage
    distance_matrix = np.loadtxt("./data/distance_matrix.txt", delimiter='\t')
    linkage_matrix = linkage(squareform(distance_matrix), method="average") # We use squareform to condense our distance matrix

    cutting_points = linkage_matrix[:, 2] # The second col of linkage_metric is all merge distances
    results = []
    result = ()
    best_result = ()
    best_silhouette = -1
    for cut in cutting_points:
        clusters = fcluster(linkage_matrix, cut, criterion = "distance")
        cluster_count = len(np.unique(clusters))

        if 1 < cluster_count < len(distance_matrix):
            
            silhouette = silhouette_score(distance_matrix, clusters, metric="precomputed")
            
            result = (
                cut,
                cluster_count,
                silhouette
            )

            if silhouette > best_silhouette:
                best_result = result
                best_silhouette = silhouette
                
            results.append(result)
            
    print(f"Best cut at distance = {best_result[0]}\nNumber of Clusters = {best_result[1]}\nSilhouette Score = {best_result[2]}")

    # Cluster Count Vs. Silhouette Score
    plt.figure()
    plt.plot([data[1] for data in results], [data[2] for data in results])
    plt.title("Cluster Count Vs. Silhouette Score")
    plt.ylabel("Silhouette Score")
    plt.xlabel("Number of Clusters")
    plt.xticks(range(0, 140, 10))
    plt.savefig('./data/silhouette_cluster_count.svg', bbox_inches='tight')

    # Tree
    plt.figure()
    dendrogram(linkage_matrix, color_threshold = 0, leaf_font_size = 3)
    plt.title("Phylogenetic Tree of 135 Bacterial Sequences")
    plt.ylabel("Distance")
    plt.xlabel("Baterica ID")
    plt.gcf().set_size_inches(8, 4.5)
    plt.savefig('./data/tree.svg', bbox_inches='tight')
