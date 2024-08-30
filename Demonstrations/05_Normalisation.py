import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import FunctionTransformer
from scipy import stats
import joblib

# Set number of cores for processing
from joblib import Parallel, delayed

# Load data
sce = joblib.load("R_objects/Caron_filtered.500.rds")

# Check data object
print(sce)

# Check number of cells per sample
print(sce['SampleName'].value_counts())

# PBMMC_1 cell UMI counts distribution before normalization
oneSamTab = sce[sce['SampleName'] == "PBMMC_1"][['SampleName', 'Barcode', 'sum']]
oneSamTab['cell_num'] = np.arange(1, len(oneSamTab) + 1)

# Make a bar plot of the UMI counts per cell
p_before_nom = sns.barplot(data=oneSamTab, x='cell_num', y='sum')
plt.xlabel('Cell Index')
plt.ylabel('Cell UMI counts')
plt.title("PBMMC_1: Before Normalization")
plt.xticks(rotation=90)
plt.show()

# Deconvolution Normalisation
np.random.seed(100)

# Get clusters - deconvolution pools are made only within clusters
clust = quickCluster(sce, n_jobs=7)

# How many clusters do we have and how big are they?
print(pd.Series(clust).value_counts())

# We have just the raw counts assay at the moment
print(sce.assay_names)

# We have no size factors for the cells at this stage
print(sce.size_factors)

# Compute pooled size factors
sce = computePooledFactors(sce, clusters=clust, min_mean=0.1, n_jobs=7)

# We still have just the raw counts assay
print(sce.assay_names)

# We now have size factors which we can extract and summarise
deconv_sf = sce.size_factors
print(deconv_sf.describe())

# Calculate library size factors
lib_sf = librarySizeFactors(sce)

# Combine the size factor results with the sample group names and make a scatter plot of the results
size_factors_df = pd.DataFrame({
    'LibrarySizeFactors': lib_sf,
    'DeconvolutionSizeFactors': deconv_sf,
    'SampleGroup': sce['SampleGroup']
})

sns.scatterplot(data=size_factors_df, x='LibrarySizeFactors', y='DeconvolutionSizeFactors', hue='SampleGroup')
plt.plot([size_factors_df['LibrarySizeFactors'].min(), size_factors_df['LibrarySizeFactors'].max()],
         [size_factors_df['LibrarySizeFactors'].min(), size_factors_df['LibrarySizeFactors'].max()], color='red')
plt.show()

# Deconvolution: scale and transform raw counts and save object
sce = logNormCounts(sce)

# PBMMC_1 cell UMI counts distribution after normalization
norm_counts = logNormCounts(sce, transform='none').assay('normcounts').sum(axis=1)

# Make a DataFrame of the total normalised counts per cell
norm_counts_df = pd.DataFrame({
    'Barcode': norm_counts.index,
    'normCounts': np.log2(norm_counts)
})

# Add the raw counts sum, cell numbers and SampleName to the DataFrame
norm_counts_df = norm_counts_df.merge(oneSamTab, on='Barcode')

# Plot the raw and normalised cell UMI counts
p_after_norm = sns.barplot(data=norm_counts_df, x='cell_num', y='normCounts')
plt.xlabel('Cell Index')
plt.ylabel('Normalized Cell UMI counts')
plt.title("PBMMC_1: After Normalization")
plt.xticks(rotation=90)
plt.show()

# Let's separate out the scaling normalisation from the log transformation
p_before_nom_nolog = sns.barplot(data=oneSamTab, x='cell_num', y=np.log2(oneSamTab['sum']))
plt.xlabel('Cell Index')
plt.ylabel('Cell UMI counts')
plt.title("Logged raw counts")
plt.xticks(rotation=90)
plt.show()

# Mean and variance for raw counts
mean = np.mean(sce['counts'], axis=0)
var = np.var(sce['counts'], axis=0)

# There is a strong linear relationship between mean and variance
plt.scatter(np.log(mean), np.log(var))
plt.plot(np.log(mean), np.log(mean), color="red")
plt.show()

# Mean and variance for scaled counts
mean_scaled = logNormCounts(sce, transform='none').assay('normcounts').mean(axis=0)
var_scaled = logNormCounts(sce, transform='none').assay('normcounts').var(axis=0)

plt.scatter(np.log(mean_scaled), np.log(var_scaled))
plt.plot(np.log(mean_scaled), np.log(mean_scaled), color="red")
plt.show()

# Mean and variance for scaled, log transformed counts
mean_norm = np.mean(sce['logcounts'], axis=0)
var_norm = np.var(sce['logcounts'], axis=0)

plt.scatter(mean_norm, var_norm)
plt.plot(mean_norm, mean_norm, color="red")
plt.show()

# Save sce object after normalisation
joblib.dump(sce, "results/caron_normalised.rds")

# Exercise 1
# Exercise: apply the deconvolution normalization on a single sample: ETV6-RUNX1_1 (aka GSM3872434).

# sctransform: Variant Stabilising Transformation
counts = sce['counts']
print(type(counts))
print(counts.shape)

# Add gene attributes
gene_attr = pd.DataFrame({
    'mean': np.mean(counts, axis=0),
    'detection_rate': np.mean(counts > 0, axis=0),
    'var': np.var(counts, axis=0)
})
gene_attr['log_mean'] = np.log10(gene_attr['mean'])
gene_attr['log_var'] = np.log10(gene_attr['var'])

print(gene_attr.shape)
print(gene_attr.head())

# Add cell attributes
cell_attr = pd.DataFrame({
    'n_umi': np.sum(counts, axis=1),
    'n_gene': np.sum(counts > 0, axis=1)
})

print(cell_attr.shape)
print(cell_attr.head())

# Over dispersion plot
sns.scatterplot(data=gene_attr, x='log_mean', y='log_var', alpha=0.3)
plt.plot([gene_attr['log_mean'].min(), gene_attr['log_mean'].max()],
         [gene_attr['log_mean'].min(), gene_attr['log_mean'].max()], color='red')
plt.show()

# Detection rate plot
x = np.linspace(-3, 2, 1000)
poisson_model = pd.DataFrame({
    'log_mean': x,
    'detection_rate': 1 - stats.poisson.pmf(0, np.power(10, x))
})

sns.scatterplot(data=gene_attr, x='log_mean', y='detection_rate', alpha=0.3)
sns.lineplot(data=poisson_model, x='log_mean', y='detection_rate', color='red')
plt.show()

# n_gene vs n_umi (Total UMI counts in a cell)
sns.scatterplot(data=cell_attr, x='n_umi', y='n_gene', alpha=0.3)
sns.kdeplot(data=cell_attr, x='n_umi', y='n_gene', levels=10, color='black', linewidths=0.3)
plt.show()

# sctransform (VST run)
np.random.seed(44)
vst_out = vst(umi=counts, latent_var=['log_umi'], return_gene_attr=True, return_cell_attr=True, verbosity=2)

# sctransform plot model parameters
plot_model_pars(vst_out, verbosity=1)

# sctransform model
print(vst_out['model_str'])

# Inspect model using two genes
gene_symbols = ['RPL10', 'HBB']
ensId = rowData(sce)[rowData(sce)['Symbol'].isin(gene_symbols)]['ID']

plot_model(x=vst_out, umi=counts, goi=ensId, plot_residual=True)

# Properties of transformed data
# Residual mean is centered around 0
sns.histplot(vst_out['gene_attr']['residual_mean'], bins=30)
plt.show()

# Properties of transformed data
# Residual variance is centered around 1
sns.histplot(vst_out['gene_attr']['residual_variance'], bins=30)
plt.axvline(x=1, color='red')
plt.xlim(0, 10)
plt.show()

# Properties of transformed data
# Mean UMI vs variance plot
sns.scatterplot(data=vst_out['gene_attr'], x='log10(gmean)', y='residual_variance', alpha=0.3)
plt.show()

# Top genes with high residual variance
top_genes = vst_out['gene_attr'].nlargest(10, 'residual_variance')
top_genes = top_genes.round(2).reset_index()
top_genes = top_genes.merge(rowData(sce)[['ID', 'Symbol']], on='ID')

# Add VST to sce
keepGenes = sce.index.isin(vst_out['y'].index)
sce = sce[keepGenes]
vstMat = vst_out['y'][sce.index].to_sparse()

sce['sctrans_norm'] = vstMat

