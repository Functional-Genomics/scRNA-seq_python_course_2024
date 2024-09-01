# Setup & Data ----

# load packages
library(scater) 
library(scran)
library(PCAtools)
library(tidyverse)

# read data
sce <- readRDS("R_objects/Caron_normalized.500.rds")
sce

# use gene symbols as rownames - making sure they are unique
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)


# Gene variance & HVGs ----

# model gene variance
gene_var <- modelGeneVar(sce)

gene_var

# plot gene variable
gene_var %>% 
  as.data.frame() %>% 
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")

# get the most highly variable genes (HGVs)
hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs)
hvgs[1:10]

# visualise their expression
plotExpression(sce, features = hvgs[1:20], point_alpha = 0.05)


# PCA ----

# this returns a new SCE object with a "PCA" reducedDim slot 
sce <- runPCA(sce, subset_row = hvgs)
sce

# extract a few rows and columns of the PCA matrix
reducedDim(sce, "PCA")[1:10, 1:5]

# extract the % variance explained by each PC
percent.var <- attr(reducedDim(sce), "percentVar")

# visualise as a "scree plot"
plot(percent.var, log = "y", xlab = "PC", ylab = "Variance explained (%)")

# visualise PCA result - PC1 vs PC2
plotReducedDim(sce, dimred = "PCA", colour_by = "SampleName")

# grid plot for multiple PCs
plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "SampleName")

# Highly customisable plots with ggcells
ggcells(sce, aes(x = PCA.1, y = PCA.2, colour = SampleName)) +
  geom_point(size = 0.5) +
  facet_wrap(~ SampleName) +
  labs(x = "PC1", y = "PC2", colour = "Sample")

# PCA diagnostics ----

# "correlation" between certain variables and our PCs
explain_pcs <- getExplanatoryPCs(sce,
                                variables = c("sum",
                                              "detected",
                                              "SampleGroup",
                                              "SampleName",
                                              "subsets_Mito_percent")
                                )

plotExplanatoryPCs(explain_pcs/100)

# "correlation" between our variables and each gene
plotExplanatoryVariables(sce,
                         variables = c(
                           "sum",
                           "detected",
                           "SampleGroup",
                           "SampleName",
                           "subsets_Mito_percent"
                         ))


# Choosing PCs ----

# elbow method
chosen_elbow <- findElbowPoint(percent.var)
chosen_elbow

# visualise the cut point on the scree plot
plot(percent.var)
abline(v=chosen_elbow, col="dodgerblue")

# "denoising" PCA - remove PCs that explain less than technical variance
sce.denoised <- denoisePCA(sce, technical = gene_var, subset.row = hvgs)

ncol(reducedDim(sce.denoised, "PCA"))


# Exercise 1: t-SNE ----

# add the t-SNE result to the reducedDim slot of the SCE object
# we name this reducedDim "TSNE_perplex50"
# we set perplexity = 50 (which is the default if we don't specify it)
# we run t-SNE based on the PCA we ran previously
# we will use the first 10 principal components
set.seed(123) # set a random seed to ensure reproducibility
sce <- runTSNE(sce,
               name = "TSNE_perplex50",
               perplexity = 50,
               dimred = "PCA",
               n_dimred = 10)

# Make a custom visualisation using ggcells
ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2,
                 colour = SampleName)) +
  geom_point()

# Part A
# Re-run the algorithm but change the random seed number.
# Do the results change dramatically between runs?
FIXME

# Part B
# Instead of colouring by SampleName, colour by expression of known cell markers
# CD79A (B cells)
# CST3 (monocytes)
# CD3D (T cells)
# HBA1 (erythrocytes)
FIXME

# Part C
# Facet these plots by SampleName to better understand where each marker is mostly expressed
FIXME

# Part D
# Explore different perplexity values (for example 5 and 500)
# Do you get tighter or looser clusters?
FIXME


# Exercise 2: UMAP ----

# Part A
# run the UMAP with 50 neighbours
set.seed(123) # set seed for reproducibility
sce <- runUMAP(sce,
               name = "UMAP_neighbors50",
               dimred = "PCA",
               FIXME)

# Part B
# visualise the resulting UMAP projection (colour cells by sample)
FIXME

# Part C
# run the UMAP with 5 and 500 neighbours and compare the results
FIXME

# Part D
# compare the UMAP projection with the t-SNE projections
# would you prefer one over the other?
FIXME
