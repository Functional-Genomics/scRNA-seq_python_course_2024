# Setup ----

# Import required libraries
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# Load the AnnData object
sce = sc.read("R_objects/Caron_clustered.PBMMCandETV6RUNX1.h5ad")

# Plot UMAP done on the batch-corrected data
sc.pl.umap(sce, color="label", legend_loc="on data", title="UMAP (batch-corrected)")

# Create pseudo-bulk samples -----

# Tabulate number of cells per label + sample combination
print(pd.crosstab(sce.obs["label"], sce.obs["SampleName"]))

# Sum counts across cells - by label and sample
def aggregate_across_cells(adata, group_by):
    grouped = adata.obs.groupby(group_by)
    aggregated = adata.X.T.dot(grouped.transform(lambda x: x == x.name).astype(float)).T
    obs = grouped.first()
    obs["ncells"] = grouped.size()
    var = adata.var.copy()
    return ad.AnnData(X=aggregated, obs=obs, var=var)

summed = aggregate_across_cells(sce, ["label", "SampleName"])

# The output is a new AnnData object with aggregated counts matrix
print(summed)

# Run DE analysis -----

# Filter our pseudo-bulk object
summed = summed[summed.obs["ncells"] >= 20]

# Perform differential analysis for each label
def pseudo_bulk_dge(adata, label, design, coef, condition):
    results = {}
    for l in adata.obs[label].unique():
        subset = adata[adata.obs[label] == l]
        design_matrix = sm.add_constant(pd.get_dummies(subset.obs[design]))
        model = sm.GLM(subset.X, design_matrix, family=sm.families.NegativeBinomial())
        result = model.fit()
        coef_index = list(design_matrix.columns).index(coef)
        lfc = result.params[coef_index]
        pvals = result.pvalues[coef_index]
        fdr = multipletests(pvals, method='fdr_bh')[1]
        results[l] = pd.DataFrame({
            'logFC': lfc,
            'PValue': pvals,
            'FDR': fdr
        }, index=adata.var_names)
    return results

de_results = pseudo_bulk_dge(summed, 'label', 'SampleGroup', 'SampleGroupPBMMC', 'SampleName')

# Extract one of the tables from the dictionary
b_c1_res = de_results['B_C1']

# Plot mean-CV relationship across genes
plt.figure(figsize=(10, 6))
plt.scatter(np.log10(b_c1_res.index.map(lambda x: summed[:, x].X.mean())), 
            np.log10(b_c1_res.index.map(lambda x: summed[:, x].X.std() / summed[:, x].X.mean())),
            alpha=0.1)
plt.xlabel('log10(mean)')
plt.ylabel('log10(CV)')
plt.title('Mean-CV relationship')
plt.show()

# Plot mean-difference plot
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for i, ax in enumerate(axes.flat):
    if i < summed.n_obs:
        sc.pl.scatter(summed, x='n_counts', y=summed.var_names[i], ax=ax, show=False)
        ax.axhline(y=0, color='salmon', linewidth=2)
    else:
        ax.axis('off')
plt.tight_layout()
plt.show()

# Plot MDS
sc.pp.normalize_total(summed, target_sum=1e6)
sc.pp.log1p(summed)
sc.tl.pca(summed)
sc.pl.pca(summed, color='SampleGroup', palette={'PBMMC': 'tomato', 'ETV6RUNX1': 'steelblue'})

# Explore DEGs ----

# Identify DEGs based on FDR threshold
def decide_tests_per_label(de_results, threshold):
    return pd.DataFrame({label: (res['FDR'] < threshold).astype(int) for label, res in de_results.items()})

is_de = decide_tests_per_label(de_results, 0.05)

# Summarize the results
print(is_de.sum())

# Filter the results table for DEGs
print(b_c1_res.sort_values('FDR').head())

# Visualize from our summed single cell object
sc.pl.violin(summed, ['HTR1F'], groupby='SampleName', rotation=45, 
             stripplot=True, jitter=0.4, size=1, palette=['tomato', 'steelblue'])

# Note: The exercise part is not included as it requires loading additional data and performing similar analyses.

