### WARNING: THIS IS A RAW CONVERSION ###

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from scipy import stats
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import scanpy.external as sce

# Import required libraries
DropletUtils = importr('DropletUtils')
scater = importr('scater')
ensembldb = importr('ensembldb')
AnnotationHub = importr('AnnotationHub')
BiocParallel = importr('BiocParallel')
tidyverse = importr('tidyverse')
ggvenn = importr('ggvenn')

# Read sample sheet
samplesheet = pd.read_csv("Data/sample_sheet.tsv", sep="\t")

# Set up parallel processing
bp_params = BiocParallel.MulticoreParam(workers=7)

# Read 10x data
sample_path = "Data/CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
sce_sing = sc.read_10x_mtx(sample_path)

# Print dimensions and first 10x10 of counts
print(sce_sing.shape)
print(sce_sing.X[:10, :10].toarray())

# Print row and column data
print(sce_sing.var)
print(sce_sing.obs)

# Count non-zero genes
print(np.sum(sce_sing.X.sum(axis=0) > 0))

# Plot genes per cell
genes_per_cell = np.sum(sce_sing.X > 0, axis=1)
sns.kdeplot(genes_per_cell)
plt.xlabel("Genes per cell")
plt.show()

# Plot mean UMIs per cell vs proportion of cells expressing the gene
mean_umis = sce_sing.X.mean(axis=0)
prop_cells = np.mean(sce_sing.X > 0, axis=0)
plt.scatter(mean_umis, prop_cells)
plt.xscale('log')
plt.xlabel("Mean UMIs per cell")
plt.ylabel("Proportion of cells expressing the gene")
plt.show()

# Plot relative expression of top 20 genes
rel_expression = sce_sing.X / sce_sing.X.sum(axis=1)[:, np.newaxis] * 100
most_expressed = rel_expression.sum(axis=0).argsort()[-20:][::-1]
plot_data = rel_expression[:, most_expressed].T
plt.boxplot(plot_data, vert=False)
plt.xlabel("% total count per cell")
plt.show()

# Read multiple samples
samples = samplesheet['Sample'].iloc[[0, 4, 6, 8]]
list_of_files = [f"Data/CellRanger_Outputs/{sample}/outs/filtered_feature_bc_matrix" for sample in samples]
adata_list = [sc.read_10x_mtx(path) for path in list_of_files]
adata = AnnData.concatenate(*adata_list, join='outer', batch_key='Sample')

# Add sample information
adata.obs = adata.obs.merge(samplesheet, on='Sample', how='left')

# Filter genes
adata = adata[:, adata.X.sum(axis=0) > 0]

# Add gene annotations
ah = AnnotationHub.AnnotationHub()
ens_hs_107 = ah.query(["Homo sapiens", "EnsDb", "107"])[[1]]
gene_annot = ensembldb.select(ens_hs_107, keys=adata.var_names, keytype="GENEID", columns=["GENEID", "SEQNAME"])
gene_annot.columns = ["ID", "Chromosome"]
adata.var = adata.var.merge(gene_annot, left_index=True, right_on="ID")

# Calculate QC metrics
adata.var['mt'] = adata.var['Chromosome'] == "MT"
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

# Identify outliers
adata.obs['outlier'] = (
    stats.zscore(np.log1p(adata.obs['n_genes_by_counts'])) < -3
) | (
    stats.zscore(np.log1p(adata.obs['total_counts'])) < -3
) | (
    stats.zscore(adata.obs['pct_counts_mt']) > 3
)

# Plot QC metrics with outliers highlighted
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             jitter=0.4, multi_panel=True, groupby='outlier')

# Filter cells
adata = adata[~adata.obs['outlier'], :]

# Final QC plot
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='SampleName')

