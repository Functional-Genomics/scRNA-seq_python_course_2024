# Import necessary libraries
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# Read single cell object
sce = sc.read("R_objects/Caron_clustered.500.h5ad")

# Check labels are set to our clusters
assert all(sce.obs['k.25_cluster.fun.leiden'] == sce.obs['label'])

# Visualize UMAP of our clusters
sc.pl.umap(sce, color='label', legend_loc='on data', title='UMAP of clusters')

# Visualize monocyte-specific marker
sc.pl.umap(sce, color='CST3', use_raw=False, title='Expression of CST3')

# Calculate pairwise marker gene statistics
sc.tl.rank_genes_groups(sce, 'label', method='wilcoxon', key_added='wilcoxon')

# Extract results for cluster 11
c11_markers = sc.get.rank_genes_groups_df(sce, group='11', key='wilcoxon').set_index('names')

# Filter markers based on rank statistics
c11_top_markers = c11_markers[c11_markers['logfoldchanges'].abs() > 0.5].sort_values('scores', ascending=False).head()

# Visualize one of the markers
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.umap(sce, color='LYZ', ax=ax1, show=False, title='Expression of LYZ')
sc.pl.violin(sce, keys='LYZ', groupby='label', ax=ax2, show=False)
plt.tight_layout()
plt.show()

# Exercise 1
sc.pl.umap(sce, color='CD3D', legend_loc='on data', title='Expression of CD3D')
sc.pl.violin(sce, keys='CD3D', groupby='label')

# Extract results for cluster 6
c6_markers = sc.get.rank_genes_groups_df(sce, group='6', key='wilcoxon').set_index('names')

# Filter the dataframe
c6_top_markers = c6_markers[c6_markers['logfoldchanges'].abs() > 0.5].sort_values('scores', ascending=False).head(10)

# Visualize expression of top genes
sc.pl.umap(sce, color=c6_top_markers.index.tolist(), ncols=3)

# Heatmaps
c11_top_genes = c11_top_markers.index.tolist()

sc.pl.heatmap(sce, var_names=c11_top_genes, groupby='label', show_gene_labels=True)

sc.pl.dotplot(sce, var_names=c11_top_genes, groupby='label')

# Adjusting log-fold change
print(c11_markers.loc['SNX10', ['logfoldchanges', 'scores']])

sc.pl.violin(sce, keys='SNX10', groupby='label')

# Recalculate markers with higher logfc threshold
sc.tl.rank_genes_groups(sce, 'label', method='wilcoxon', key_added='wilcoxon_lfc2', logfc_threshold=2)

c11_markers_lfc2 = sc.get.rank_genes_groups_df(sce, group='11', key='wilcoxon_lfc2').set_index('names')
print(c11_markers_lfc2.loc['SNX10', 'scores'])

# Annotation labels
known_genes = ['HBA1', 'CST3', 'CD3E', 'NKG7', 'CD79A', 'MS4A1']

sc.pl.violin(sce, keys=known_genes, groupby='label', rotation=90)

sc.pl.heatmap(sce, var_names=known_genes, groupby='label', show_gene_labels=True)

# Re-label the cells
new_labels = {
    '1': 'B (c1)', '2': 'B (c2)', '3': 'B (c3)', '4': 'B (c4)',
    '5': 'CD20+ B (c5)', '6': 'T (c6)', '7': 'NK T (c7)',
    '8': 'Erythrocytes (c8)', '9': 'Erythrocytes (c9)',
    '10': 'Erythrocytes (c10)', '11': 'Monocytes (c11)', '12': 'B (c12)'
}

sce.obs['new_label'] = sce.obs['label'].map(new_labels)

# Visualize UMAP with new labels
sc.pl.umap(sce, color='new_label', legend_loc='on data', title='UMAP with new labels')

