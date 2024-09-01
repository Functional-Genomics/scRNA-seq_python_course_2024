# Import necessary libraries
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import LabelEncoder
from sklearn.manifold import TSNE
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csr_matrix
import umap

# Set working directory
import os
os.chdir("~/Course_Materials")

# Load data
sce = sc.read_rds("R_objects/Caron_batch_corrected.500.rds")

# Check samples
sample_counts = sce.obs['SampleName'].value_counts()

###############################
# Clustering with the walktrap method
#####################################

# Cluster cells using default clustering parameters and the batch-corrected reduced dimensions
clustering1 = sc.tl.leiden(sce, resolution=1.0)

# Table of clusters
cluster_counts = clustering1.value_counts()

# Plot tSNE with default clusters
sce.obs['Clusters1'] = clustering1
sc.tl.tsne(sce, use_rep='corrected')
sc.pl.tsne(sce, color='Clusters1')

# Generate an alternative clustering using Shared Nearest Neighbours with k = 15
sce.obs['walktrap15'] = sc.tl.leiden(sce, resolution=1.0, key_added='walktrap15')

# A heatmap showing the (logged) numbers of cells in each cluster for each dataset with walktrap k = 15
w15_table = np.log(pd.crosstab(sce.obs['walktrap15'], sce.obs['SampleName']) + 1)
sns.heatmap(w15_table, cmap='viridis')

# Plot tSNE walktrap k15
sc.pl.tsne(sce, color='walktrap15')

# Plot tSNE walktrap k15 by sample group
sc.pl.tsne(sce, color='walktrap15', add_outline=True, group='SampleGroup')

###############################
# Clustering with the Louvain method
#####################################

# Cluster cells - louvain k15
sce.obs['louvain15'] = sc.tl.leiden(sce, resolution=1.0, key_added='louvain15')

# Plot tSNE louvain k15
sc.pl.tsne(sce, color='louvain15')

# Plot tSNE louvain k15 by sample group
sc.pl.tsne(sce, color='louvain15', add_outline=True, group='SampleGroup')

####################################
# Assessing cluster behaviour
###################################

# Calculate silhouette widths for the Leiden, k=20 clustering
silhouette_vals = silhouette_samples(sce.obsm['X_corrected'], sce.obs['leiden20'])

# Examine the silhouette data
silhouette_df = pd.DataFrame(silhouette_vals, columns=['silhouette_width'])
silhouette_df['cluster'] = sce.obs['leiden20']

# Silhouette width plot with beeswarm
plt.figure(figsize=(10, 6))
sns.swarmplot(x='cluster', y='silhouette_width', data=silhouette_df, color='blue', alpha=0.6)
plt.title('Silhouette Widths for Clusters')
plt.show()

# Silhouette width grid
def plot_silhouette_grid(silDat):
    silDat['closestCluster'] = np.where(silDat['silhouette_width'] > 0, silDat['cluster'], 'other')
    grid_data = silDat.groupby(['cluster', 'closestCluster']).size().reset_index(name='olap')
    grid_data['total'] = grid_data.groupby('cluster')['olap'].transform('sum')
    grid_data['proportion'] = grid_data['olap'] / grid_data['total']
    grid_data['proportion'] = np.where(grid_data['cluster'] == grid_data['closestCluster'], grid_data['proportion'], -grid_data['proportion'])
    
    pivot_table = grid_data.pivot('cluster', 'closestCluster', 'proportion')
    sns.heatmap(pivot_table, cmap='RdYlGn', center=0)
    plt.title('Silhouette Width Grid')
    plt.show()

# Save Leiden 20 sil grid plot
plot_silhouette_grid(silhouette_df)

# Plot silhouette width for walktrap k=15 clustering
silhouette_vals_wp = silhouette_samples(sce.obsm['X_corrected'], sce.obs['walktrap15'])
silhouette_df_wp = pd.DataFrame(silhouette_vals_wp, columns=['silhouette_width'])
silhouette_df_wp['cluster'] = sce.obs['walktrap15']
plot_silhouette_grid(silhouette_df_wp)

# Plot silhouette width for louvain k=15 clustering
silhouette_vals_lp = silhouette_samples(sce.obsm['X_corrected'], sce.obs['louvain15'])
silhouette_df_lp = pd.DataFrame(silhouette_vals_lp, columns=['silhouette_width'])
silhouette_df_lp['cluster'] = sce.obs['louvain15']
plot_silhouette_grid(silhouette_df_lp)

# Look at the concordance of different clusterings using the Jaccard index
jacc_mat = pd.crosstab(sce.obs['louvain15'], sce.obs['walktrap15'])
sns.heatmap(jacc_mat, cmap='viridis')

# Use clusterSweep to examine a range of clustering parameters
# clusterSweep with k = 5, 10, 15, 20, 25 and walktrap method
# Note: This part requires a custom implementation of clusterSweep in Python

# Make a data frame of clusterSweep results
# Note: This part requires a custom implementation of clusterSweep in Python

# Generate line plots of the cluster number and mean silhouette width for different values of k
# Note: This part requires a custom implementation of clusterSweep in Python

# Add clusterSweep output to sce
# Note: This part requires a custom implementation of clusterSweep in Python

# Set labels to our favourite clustering
sce.obs['label'] = sce.obs['k.25_cluster.fun.leiden']

# Plot tSNE leiden k25
sc.pl.tsne(sce, color='label')

# Switch the rownames in the SCE object to be gene symbols
sce.var_names = pd.Series(sce.var['Symbol']).unique()

# Plot B cell marker expression on TSNE
sc.pl.tsne(sce, color='CD79A')

# Plot expression B cell marker as a violin plot
sns.violinplot(x='label', y='CD79A', data=sce.obs)

# Plot tSNE monocyte markers
sc.pl.tsne(sce, color='LYZ')

# Plot expression monocyte markers
sns.violinplot(x='label', y='LYZ', data=sce.obs)

