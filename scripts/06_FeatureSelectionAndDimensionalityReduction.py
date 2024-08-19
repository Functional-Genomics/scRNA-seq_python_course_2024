# Setup & Data
# Import required libraries
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read data
sce = sc.read("R_objects/Caron_normalized.500.h5ad")
print(sce)

# Use gene symbols as index - making sure they are unique
sce.var_names_make_unique()

# Gene variance & HVGs

# Model gene variance
sc.pp.highly_variable_genes(sce, n_top_genes=int(0.1 * sce.n_vars), flavor="seurat_v3")

# Plot gene variable
plt.figure(figsize=(10, 6))
sc.pl.highly_variable_genes(sce)

# Get the most highly variable genes (HVGs)
hvgs = sce.var.index[sce.var.highly_variable]
print(f"Number of HVGs: {len(hvgs)}")
print(f"First 10 HVGs: {hvgs[:10].tolist()}")

# Visualise their expression
sc.pl.violin(sce, hvgs[:20], jitter=0.4, alpha=0.05, multi_panel=True)

# PCA

# Run PCA
sc.tl.pca(sce, use_highly_variable=True)

# Extract a few rows and columns of the PCA matrix
print(sce.obsm["X_pca"][:10, :5])

# Extract the % variance explained by each PC
percent_var = sce.uns["pca"]["variance_ratio"] * 100

# Visualise as a "scree plot"
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(percent_var) + 1), percent_var)
plt.yscale("log")
plt.xlabel("PC")
plt.ylabel("Variance explained (%)")
plt.show()

# Visualise PCA result - PC1 vs PC2
sc.pl.pca(sce, color="SampleName")

# Grid plot for multiple PCs
sc.pl.pca_variance_ratio(sce, n_pcs=50)

# PCA diagnostics

# "Correlation" between certain variables and our PCs
sc.pl.pca_overview(
    sce, keys=["sum", "n_genes_by_counts", "SampleGroup", "SampleName", "pct_counts_mt"]
)

# Choosing PCs

# Elbow method
from kneed import KneeLocator

kn = KneeLocator(
    range(1, len(percent_var) + 1), percent_var, curve="convex", direction="decreasing"
)
chosen_elbow = kn.elbow

print(f"Chosen elbow point: {chosen_elbow}")

# Visualise the cut point on the scree plot
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(percent_var) + 1), percent_var)
plt.axvline(x=chosen_elbow, color="dodgerblue")
plt.xlabel("PC")
plt.ylabel("Variance explained (%)")
plt.show()

# Exercise 1: t-SNE

# Add the t-SNE result
sc.tl.tsne(sce, n_pcs=10, perplexity=50, random_state=123)

# Visualize t-SNE
sc.pl.tsne(sce, color="SampleName")

# Part A: Re-run with different random seed
sc.tl.tsne(sce, n_pcs=10, perplexity=50, random_state=456, key_added="tsne_seed456")
sc.pl.tsne(sce, color="SampleName", use_rep="X_tsne_seed456")

# Part B: Color by expression of known cell markers
markers = ["CD79A", "CST3", "CD3D", "HBA1"]
sc.pl.tsne(sce, color=markers)

# Part C: Facet plots by SampleName
for marker in markers:
    sc.pl.tsne(sce, color=marker, facet_wrap="SampleName")

# Part D: Explore different perplexity values
for perplexity in [5, 500]:
    sc.tl.tsne(
        sce,
        n_pcs=10,
        perplexity=perplexity,
        random_state=123,
        key_added=f"tsne_perplex{perplexity}",
    )
    sc.pl.tsne(sce, color="SampleName", use_rep=f"X_tsne_perplex{perplexity}")

# Exercise 2: UMAP

# Part A: Run UMAP with 50 neighbors
sc.pp.neighbors(sce, n_neighbors=50, n_pcs=10)
sc.tl.umap(sce)

# Part B: Visualize UMAP projection
sc.pl.umap(sce, color="SampleName")

# Part C: Run UMAP with 5 and 500 neighbors
for n_neighbors in [5, 500]:
    sc.pp.neighbors(
        sce, n_neighbors=n_neighbors, n_pcs=10, key_added=f"neighbors_{n_neighbors}"
    )
    sc.tl.umap(
        sce, neighbors_key=f"neighbors_{n_neighbors}", key_added=f"umap_{n_neighbors}"
    )
    sc.pl.umap(sce, color="SampleName", use_rep=f"X_umap_{n_neighbors}")

# Part D: Compare UMAP and t-SNE projections
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
sc.pl.umap(sce, color="SampleName", ax=ax1, show=False)
ax1.set_title("UMAP")
sc.pl.tsne(sce, color="SampleName", ax=ax2, show=False)
ax2.set_title("t-SNE")
plt.tight_layout()
plt.show()
