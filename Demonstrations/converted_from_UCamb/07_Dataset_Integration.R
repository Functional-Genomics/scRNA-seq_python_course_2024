############################################################################
# load_packages 
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(pheatmap)
library(magrittr)
############################################################################

############################################################################
# set correct working directory
setwd("~/Course_Materials")
############################################################################

############################################################################
# load and explore data
sce_rep1 <- readRDS("R_objects/PBMMC_1a_dimRed.rds")
sce_rep2 <- readRDS("R_objects/PBMMC_1b_dimRed.rds")

dim(sce_rep1)
dim(sce_rep2)
# add batch info 
colData(sce_rep1)$batch <- "1"
colData(sce_rep2)$batch <- "2"
############################################################################

############################################################################
# fit mean variance model each batch separately
gene_var_rep1 <- modelGeneVar(sce_rep1)
gene_var_rep2 <- modelGeneVar(sce_rep2)
############################################################################

############################################################################
# prepare data: get common genes between two batches and subset 
common_genes <- intersect(rownames(sce_rep1), rownames(sce_rep2))
length(common_genes)
# Subset the SCE object
sce_rep1 <- sce_rep1[common_genes, ]
dim(sce_rep1)
sce_rep2 <- sce_rep2[common_genes, ]
dim(sce_rep2)

# Subset the mean-variance results
gene_var_rep1 <- gene_var_rep1[common_genes, ]
gene_var_rep2 <- gene_var_rep2[common_genes, ]

dim(gene_var_rep1)
dim(gene_var_rep2)
############################################################################

############################################################################
# Re-scale the data between two batches
rescaled_sces <- multiBatchNorm(sce_rep1, sce_rep2)

class(rescaled_sces)
length(rescaled_sces)
sce <- cbind(rescaled_sces[[1]], 
             rescaled_sces[[2]])
############################################################################

############################################################################
# combine highly variable genes from two batches
gene_var_combined <- combineVar(gene_var_rep1, gene_var_rep2)

# Select highly variable genes 
hvgs <- gene_var_combined$bio > 0
############################################################################

############################################################################
# run PCA on Re-scale 
sce <- runPCA(sce, subset_row = hvgs)

# cluster cells
sce$cluster_uncorrected <- clusterCells(sce, use.dimred = "PCA")

# generate TSNE space
sce <- runTSNE(sce, dimred = "PCA", name = "TSNE_uncorrected")

# Visualize uncorrected data using tSNE 
plotReducedDim(sce, dimred = "TSNE_uncorrected",
               colour_by = "batch",
               text_by = "cluster_uncorrected")


# assess cluster occupancy data as a table. 
table(Cluster = sce$cluster_uncorrected, Batch = sce$batch)
############################################################################

############################################################################
# Correcting the data using MNN

mnn_corrected <- fastMNN(sce, 
                         batch = sce$batch,
                         d = 50,
                         k = 20, 
                         subset.row = hvgs)

# explore mnn_corrected object
mnn_corrected
assay(mnn_corrected, "reconstructed")
reducedDim(mnn_corrected, "corrected")

# extract MNN corrected data 
reducedDim(sce, "corrected") <- reducedDim(mnn_corrected, "corrected")
############################################################################

############################################################################

# cluster cells using corrected PCA
sce$cluster_corrected <- clusterCells(sce, use.dimred = "corrected")

# run TSNE using corrected PCA space
sce <- runTSNE(sce, dimred = "corrected", name = "TSNE_corrected")
# Visualise corrected data using tSNE 
plotReducedDim(sce, 
               dimred = "TSNE_corrected", 
               colour_by = "batch", 
               text_by = "cluster_corrected")
############################################################################

############################################################################
# Visualise the proportion of each batch in different clusters before and after 
# correction
# uncorrected bar-plot
data.frame(Cluster = sce$cluster_uncorrected, Batch = sce$batch) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Batch), position = "fill") +
    labs(title = "MNN-uncorrected data") +
  scale_y_continuous(labels = scales::percent)

# corrected bar-plot
data.frame(Cluster = sce$cluster_corrected, Batch = sce$batch) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Batch), position = "fill") +
    labs(title = "MNN-corrected data") +
  scale_y_continuous(labels = scales::percent)
############################################################################

############################################################################
# Using the wrapper function quickCorrect 
?quickCorrect
sce_quick_mnn <- quickCorrect(sce_rep1, sce_rep2)$corrected

sce_quick_mnn$batch <- factor(sce_quick_mnn$batch)
sce_quick_mnn %>%
  runTSNE(dimred = "corrected") %>%
  plotTSNE(colour_by = "batch")
############################################################################

############################################################################
## EXERCISE 1 ## Multiple Batches
# load the data
sce_all <- readRDS("R_objects/Caron_dimRed.500.rds")

# tabulate the number of cells per sample
table(sce_all$SampleName)


# obtain a batch-corrected SCE object
sce_all_corrected <- quickCorrect(sce_all, batch = sce_all$SampleName)$corrected

# add the corrected matrix to the original object - to keep it all together
reducedDim(sce_all, "corrected") <- reducedDim(sce_all_corrected, "corrected")

#  add a tSNE using the corrected data
set.seed(323)
sce_all <- runTSNE(sce_all, 
                   dimred = "corrected",
                   name = "TSNE_corrected")

# visualise both corrected and uncorrected
plotReducedDim(sce_all, dimred = "TSNE", colour_by = "SampleName")

plotReducedDim(sce_all, dimred = "TSNE_corrected", colour_by = "SampleName")

############################################################################

############################################################################
# Applying a merge order during correction

merge_order <- list(
          list("ETV6-RUNX1_1", "ETV6-RUNX1_2", "ETV6-RUNX1_3", "ETV6-RUNX1_4"),
          list("HHD_1", "HHD_2"),
          list("PBMMC_1", "PBMMC_2", "PBMMC_3"),
          list("PRE-T_1", "PRE-T_2")
                   )

sce_all_corrected <- quickCorrect(sce_all,
                                  batch = sce_all$SampleName,
                                  PARAM = FastMnnParam(merge.order = merge_order)
                                  )


sce_all_corrected <- sce_all_corrected$corrected

# add the corrected matrix to the original object - to keep it all together
reducedDim(sce_all, "corrected_mo") <- reducedDim(sce_all_corrected, "corrected")

#  add a tSNE using the merge order corrected data
set.seed(323)
sce_all <- runTSNE(sce_all, 
                   dimred = "corrected_mo",
                   name = "TSNE_corrected_mo")


# Uncorrected tSNE 
plotReducedDim(sce_all, dimred = "TSNE", colour_by = "SampleName")


# Corrected tSNE
plotReducedDim(sce_all, dimred = "TSNE_corrected", colour_by = "SampleName")


# Merge order corrected tSNE 
plotReducedDim(sce_all, dimred = "TSNE_corrected_mo", colour_by = "SampleName")

############################################################################

############################################################################
## Correction Diagnostics 

# Mixing between batches - bar plots

# using uncorrected data
sce_all$cluster_uncorrected <- clusterCells(sce_all, use.dimred = "PCA")

data.frame(Cluster = sce_all$cluster_uncorrected, Sample = sce_all$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "MNN-Uncorrected data") +
  scale_y_continuous(labels = scales::percent)

# using corrected data
sce_all$cluster_corrected <- clusterCells(sce_all, use.dimred = "corrected")

data.frame(Cluster = sce_all$cluster_corrected, Sample = sce_all$SampleName) %>%
  ggplot(aes(x = Cluster)) +
    geom_bar(aes(fill = Sample), position = "fill") +
    labs(title = "MNN-corrected data") +
  scale_y_continuous(labels = scales::percent)

############################################################################

############################################################################
# Mixing between batches - cluster variance

cluster_var <- clusterAbundanceVar(sce_all$cluster_corrected, 
                                   batch = sce_all$SampleName)
cluster_var

batch_per_cluster <- table(Cluster = sce_all$cluster_corrected, 
                           Batch = sce_all$SampleName)
batch_per_cluster
batch_per_cluster[order(cluster_var, decreasing = TRUE), ]
############################################################################

############################################################################
## Preserving Biological Heterogeneity
# Nesting of clusters before and after correction
# using two sample (batch data)
# sce = sce_rep1 + sce_rep2
original_clusters <- colLabels(sce_rep1)
corrected_clusters_rep1 <- sce[,colnames(sce_rep1)]$cluster_corrected
tab <- nestedClusters(ref=paste("before", original_clusters),
                      alt=paste("after", corrected_clusters_rep1))

# cluster heatmap rep1 
pheatmap(tab$proportions, 
         cluster_row=FALSE, 
         cluster_col=FALSE,
         main="Sample 1 comparison")


# cluster heatmap rep2 
original_clusters <- colLabels(sce_rep2)
corrected_clusters_rep2 <- sce[,colnames(sce_rep2)]$cluster_corrected
tab <- nestedClusters(ref=paste("before", original_clusters),
                      alt=paste("after", corrected_clusters_rep2))

pheatmap(tab$proportions,
         cluster_row=FALSE,
         cluster_col=FALSE,
         main="Sample 2 comparison")

############################################################################

############################################################################
# Adjusted Rand Index - overall

pairwiseRand(corrected_clusters_rep1, colLabels(sce_rep1), mode="index")
pairwiseRand(corrected_clusters_rep1, colLabels(sce_rep1), mode="index")

pairwiseRand(corrected_clusters_rep2, colLabels(sce_rep2), mode="index")


# Adjusted Rand Index - cluster-wise

tab <- pairwiseRand(colLabels(sce_rep1), corrected_clusters_rep1)

pheatmap(tab, 
         cluster_row=FALSE, 
         cluster_col=FALSE,
         main="Sample 1 probabilities")

tab <- pairwiseRand(colLabels(sce_rep2), corrected_clusters_rep2)

pheatmap(tab, 
         cluster_row=FALSE, 
         cluster_col=FALSE,
         main="Sample 2 probabilities")
############################################################################

############################################################################
# fastMNN - lost variance
metadata(sce_quick_mnn)$merge.info$lost.var
############################################################################
