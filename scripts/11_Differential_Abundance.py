# Setup -----

# Import required libraries
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
import networkx as nx

# Load the AnnData object (equivalent to SCE in R)
adata = sc.read("R_objects/Caron_clustered.PBMMCandETV6RUNX1.h5ad")

# Check the contents of the object
print(adata)

# Plot UMAP done on the batch-corrected data
sc.pl.umap(adata, color="label", legend_loc="on data", title="UMAP (batch-corrected)")

# Differential abundance with Milo-like approach ----

# Create a Milo-like object (we'll use AnnData as a base)
milo = adata.copy()

# Build KNN graph ----

# Add KNN graph to Milo object
sc.pp.neighbors(milo, n_neighbors=60, n_pcs=50, use_rep="X_corrected")

# Sample index cells to define neighbourhoods
n_samples = int(0.1 * milo.n_obs)
sampled_indices = np.random.choice(milo.n_obs, size=n_samples, replace=False)
milo.uns['nhoods'] = milo.obsp['distances'][sampled_indices].tocsr()

# Distribution of neighbourhood sizes
plt.figure(figsize=(10, 5))
plt.hist(milo.uns['nhoods'].sum(axis=1), bins=50)
plt.axvline(x=100, color='salmon')
plt.title("Distribution of neighbourhood sizes")
plt.xlabel("Neighbourhood size")
plt.ylabel("Frequency")
plt.show()

# Count cells in each neighbourhood
milo.uns['nhood_counts'] = pd.DataFrame(
    milo.uns['nhoods'].dot(pd.get_dummies(milo.obs['SampleName'])),
    index=sampled_indices
)

# Run DA analysis ----

# Calculate distances between neighbourhoods - for p-value correction
milo.uns['nhood_distances'] = milo.uns['nhoods'].dot(milo.obsm['X_corrected'])

# Define a table for our model design
sample_info = milo.obs[['SampleName', 'SampleGroup']].drop_duplicates().set_index('SampleName')

# Run DA test (simplified version, you might need to implement a more sophisticated test)
from statsmodels.formula.api import ols

da_results = []
for i in range(milo.uns['nhood_counts'].shape[0]):
    counts = milo.uns['nhood_counts'].iloc[i]
    df = pd.DataFrame({'counts': counts, 'group': sample_info.loc[counts.index, 'SampleGroup']})
    model = ols('counts ~ group', data=df).fit()
    da_results.append({
        'logFC': model.params['group[T.ETV6-RUNX1]'],
        'PValue': model.pvalues['group[T.ETV6-RUNX1]']
    })

da_results = pd.DataFrame(da_results)
da_results['SpatialFDR'] = multipletests(da_results['PValue'], method='fdr_bh')[1]

# Visualizations ----

# P-value histogram
plt.figure(figsize=(10, 5))
plt.hist(da_results['PValue'], bins=50)
plt.title("P-value histogram")
plt.xlabel("P-value")
plt.ylabel("Frequency")
plt.show()

# Volcano plot
plt.figure(figsize=(10, 5))
plt.scatter(da_results['logFC'], -np.log10(da_results['SpatialFDR']), c=da_results['SpatialFDR'] < 0.1)
plt.axhline(y=1, color='r', linestyle='--')
plt.title("Volcano plot")
plt.xlabel("logFC")
plt.ylabel("-log10(SpatialFDR)")
plt.show()

# Build neighbourhood graph embedding and plot
G = nx.from_scipy_sparse_matrix(milo.uns['nhoods'])
pos = nx.spring_layout(G)

plt.figure(figsize=(15, 5))
plt.subplot(121)
sc.pl.umap(milo, color="label", title="UMAP with cell labels", show=False)
plt.subplot(122)
nx.draw(G, pos, node_color=da_results['logFC'], node_size=20, cmap='coolwarm')
plt.title("Neighbourhood graph")
plt.tight_layout()
plt.show()

# Annotate NHoods ----

# Annotate neighbourhood DA results with cell labels
nhood_labels = milo.uns['nhoods'].dot(pd.get_dummies(milo.obs['label']))
da_results['label'] = nhood_labels.idxmax(axis=1)
da_results['label_fraction'] = nhood_labels.max(axis=1) / nhood_labels.sum(axis=1)

# Histogram of fraction of cells in the neighbourhood with the same label
plt.figure(figsize=(10, 5))
plt.hist(da_results['label_fraction'], bins=50)
plt.title("Fraction of cells with the same label in neighbourhoods")
plt.xlabel("Label fraction")
plt.ylabel("Frequency")
plt.show()

# Add "mixed" label to neighbourhoods with less 70% consistency
da_results['label'] = np.where(da_results['label_fraction'] < 0.7, "Mixed", da_results['label'])

# Distribution of logFC across neighbourhood labels
plt.figure(figsize=(12, 6))
sns.boxplot(x='label', y='logFC', data=da_results)
plt.xticks(rotation=45, ha='right')
plt.title("Distribution of logFC across neighbourhood labels")
plt.tight_layout()
plt.show()

