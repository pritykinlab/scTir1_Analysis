import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
import seaborn as sns
def KMeans_Cluster(df_for_clustering, n_clusters=9):
    # Fit KMeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init='auto')
    kmeans.fit(df_for_clustering)

    # Original cluster labels
    original_labels = kmeans.labels_

    # Compute mean across rows for each cluster
    cluster_means = pd.DataFrame(df_for_clustering).assign(cluster=original_labels).groupby('cluster').mean().mean(axis=1)

    # Sort clusters from high to low mean
    sorted_clusters = cluster_means.sort_values(ascending=False).index

    # Create a mapping from original cluster to reordered cluster
    cluster_mapping = {old: new for new, old in enumerate(sorted_clusters)}

    # Apply mapping to cluster labels
    reordered_labels = pd.Series(original_labels).map(cluster_mapping)

    order = np.argsort(reordered_labels)


    # Generate a palette
    palette = sns.color_palette('tab20', n_colors=n_clusters)

    # Map reordered labels to colors
    row_colors = reordered_labels.map(lambda x: palette[x])
    row_colors.index = df_for_clustering.index

    clustered_genes = pd.Series(reordered_labels.values, index=df_for_clustering.index).sort_values()
    return order, row_colors, clustered_genes