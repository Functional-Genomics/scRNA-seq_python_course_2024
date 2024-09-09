import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# set correct working directory
import os
os.chdir("~/Course_Materials")

# load and explore data
sce_rep1 = sc.read("R_objects/PBMMC_1a_dimRed.rds")
sce_rep2 = sc.read("R_objects/PBMMC_1b_dimRed.rds")

print(sce_rep1.shape)
print(sce_rep2.shape)

# add batch info
sce_rep1.obs['batch'] = "1"
sce_rep2.obs['batch'] = "2"

# fit mean variance model each batch separately
gene_var_rep1 = sc.pp.highly_variable_genes(sce_rep1, flavor='seurat', n_top_genes=None, inplace=False)
gene_var_rep2 = sc.pp.highly_variable_genes(sce_rep2, flavor='seurat', n_top_genes=None, inplace=False)

# prepare data: get common genes between two batches and subset
common_genes = np.intersect1d(sce_rep1.var_names, sce_rep2.var_names)
print(len(common_genes))

# Subset the SCE object
sce_rep1 = sce_rep1[:, common_genes]
print(sce_rep1.shape)
sce_rep2 = sce_rep2[:, common_genes]
print(sce_rep2.shape)

# Subset the mean-variance results
gene_var_rep1 = gene_var_rep1[gene_var_rep1['highly_variable'], :]
gene_var_rep2 = gene_var_rep2[gene_var_rep2['highly_variable'], :]

print(gene_var_rep1.shape)
print(gene_var_rep2.shape)

# Re-scale the data between two batches
rescaled_sces = [sce_rep1, sce_rep2]  # Placeholder for multiBatchNorm equivalent

# Combine highly variable genes from two batches
gene_var_combined = pd.concat([gene_var_rep1, gene_var_rep2], axis=0)

# Select highly variable genes
hvgs = gene_var_combined['highly_variable'] > 0

# run PCA on Re-scale
sc.tl.pca(sce_rep1, use_highly_variable=True)

# cluster cells
sce_rep1.obs['cluster_uncorrected'] = sc.tl.leiden(sce_rep1)

# generate TSNE space
sc.tl.tsne(sce_rep1, use_rep='X_pca', random_state=0)

# Visualize uncorrected data using tSNE
sc.pl.tsne(sce_rep1, color='batch', title='Uncorrected Data')

# assess cluster occupancy data as a table
cluster_table = pd.crosstab(sce_rep1.obs['cluster_uncorrected'], sce_rep1.obs['batch'])
print(cluster_table)

# Correcting the data using MNN
mnn_corrected = None  # Placeholder for fastMNN equivalent

# explore mnn_corrected object
print(mnn_corrected)

# extract MNN corrected data
sce_rep1.obsm['X_corrected'] = mnn_corrected

# cluster cells using corrected PCA
sce_rep1.obs['cluster_corrected'] = sc.tl.leiden(sce_rep1, use_rep='X_corrected')

# run TSNE using corrected PCA space
sc.tl.tsne(sce_rep1, use_rep='X_corrected', random_state=0)

# Visualise corrected data using tSNE
sc.pl.tsne(sce_rep1, color='batch', title='Corrected Data')

# Visualise the proportion of each batch in different clusters before and after correction
uncorrected_df = pd.DataFrame({'Cluster': sce_rep1.obs['cluster_uncorrected'], 'Batch': sce_rep1.obs['batch']})
sns.countplot(data=uncorrected_df, x='Cluster', hue='Batch', position='fill')
plt.title("MNN-uncorrected data")
plt.ylabel("Proportion")
plt.show()

corrected_df = pd.DataFrame({'Cluster': sce_rep1.obs['cluster_corrected'], 'Batch': sce_rep1.obs['batch']})
sns.countplot(data=corrected_df, x='Cluster', hue='Batch', position='fill')
plt.title("MNN-corrected data")
plt.ylabel("Proportion")
plt.show()

# Using the wrapper function quickCorrect
sce_quick_mnn = None  # Placeholder for quickCorrect equivalent

# Add the corrected matrix to the original object
sce_all = sc.read("R_objects/Caron_dimRed.500.rds")

# tabulate the number of cells per sample
sample_counts = sce_all.obs['SampleName'].value_counts()
print(sample_counts)

# obtain a batch-corrected SCE object
sce_all_corrected = None  # Placeholder for quickCorrect equivalent

# add a tSNE using the corrected data
sce_all.obsm['X_corrected'] = sce_all_corrected

# add a tSNE using the corrected data
sc.tl.tsne(sce_all, use_rep='X_corrected', random_state=323)

# visualise both corrected and uncorrected
sc.pl.tsne(sce_all, color='SampleName', title='Uncorrected Data')
sc.pl.tsne(sce_all, color='SampleName', title='Corrected Data')

# Applying a merge order during correction
merge_order = [
    ["ETV6-RUNX1_1", "ETV6-RUNX1_2", "ETV6-RUNX1_3", "ETV6-RUNX1_4"],
    ["HHD_1", "HHD_2"],
    ["PBMMC_1", "PBMMC_2", "PBMMC_3"],
    ["PRE-T_1", "PRE-T_2"]
]

sce_all_corrected = None  # Placeholder for quickCorrect with merge order equivalent

# add the corrected matrix to the original object
sce_all.obsm['X_corrected_mo'] = sce_all_corrected

# add a tSNE using the merge order corrected data
sc.tl.tsne(sce_all, use_rep='X_corrected_mo', random_state=323)

# Uncorrected tSNE
sc.pl.tsne(sce_all, color='SampleName', title='Uncorrected Data')

# Corrected tSNE
sc.pl.tsne(sce_all, color='SampleName', title='Corrected Data')

# Merge order corrected tSNE
sc.pl.tsne(sce_all, color='SampleName', title='Merge Order Corrected Data')

# Mixing between batches - bar plots
sce_all.obs['cluster_uncorrected'] = sc.tl.leiden(sce_all)

uncorrected_df = pd.DataFrame({'Cluster': sce_all.obs['cluster_uncorrected'], 'Sample': sce_all.obs['SampleName']})
sns.countplot(data=uncorrected_df, x='Cluster', hue='Sample', position='fill')
plt.title("MNN-Uncorrected data")
plt.ylabel("Proportion")
plt.show()

sce_all.obs['cluster_corrected'] = sc.tl.leiden(sce_all, use_rep='X_corrected')

corrected_df = pd.DataFrame({'Cluster': sce_all.obs['cluster_corrected'], 'Sample': sce_all.obs['SampleName']})
sns.countplot(data=corrected_df, x='Cluster', hue='Sample', position='fill')
plt.title("MNN-corrected data")
plt.ylabel("Proportion")
plt.show()

# Mixing between batches - cluster variance
cluster_var = None  # Placeholder for clusterAbundanceVar equivalent
print(cluster_var)

batch_per_cluster = pd.crosstab(sce_all.obs['cluster_corrected'], sce_all.obs['SampleName'])
print(batch_per_cluster)

# Adjusted Rand Index - overall
ari_rep1 = None  # Placeholder for pairwiseRand equivalent
ari_rep2 = None  # Placeholder for pairwiseRand equivalent

# Adjusted Rand Index - cluster-wise
tab_rep1 = None  # Placeholder for pairwiseRand equivalent
tab_rep2 = None  # Placeholder for pairwiseRand equivalent

# fastMNN - lost variance
lost_var = None  # Placeholder for lost variance equivalent
print(lost_var)

