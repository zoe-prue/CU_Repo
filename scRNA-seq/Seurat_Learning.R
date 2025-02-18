# Learning About Seurat
# 03/28/2024
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)

# ways to enhance speed and performance of Seurat package
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# other packages to enhance seurat experience
# Signac: analysis of single-cell chromatin data
# SeuratData: automatically load datasets pre-packaged as Seurat objects
# Azimuth: local annotation of scRNA-seq and scATAC-seq queries across multiple organs and tissues
# SeuratWrappers: enables use of additional integration and differential expression methods
# # Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

# Getting Started with Seurat
# dataset of 2,700 Peripheral Blood Mononuclear Cells (PBMCs) made publicly available by 10X Genomics
# major components of a standard unsupervised clustering workflow including QC and data filtration, 
# calculation of high-variance genes, dimensional reduction, 
# graph-based clustering, and the identification of cluster markers\

# Read10X() function reads in the output of the cellranger pipeline from 10X, 
# returning a unique molecular identified (UMI) count matrix
# matrix represents the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
# more recent versions of cellranger output using the h5 file format, 
# which can be read in using the Read10X_h5() function in Seurat.

library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Desktop/CU_coding/scRNA-seq/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


# Standard pre-processing workflow
# selection and filtration of cells based on QC metrics, data normalization and scaling, 
# and the detection of highly variable features

# QC and selecting cells for further analysis
# qcs commonly used by the community:
# number of unique genes detected in each cell (nFeature_RNA)
  # low quality cells or empty droplets have very few genes
  # Cell doublets or multiplets may exhibit a high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes) (nCount_RNA)
# The percentage of reads that map to the mitochondrial genome
  # Low-quality / dying cells have high mitochondrial contamination
  # mitochondrial qc metrics calculated with PercentageFeatureSet() function,
  # which calcs the % of counts originating from a set of features (this case genes starting with MT-, or mito genes)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# final qc verification
# setting the paramters for QC variables
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Normalizing the data
# employ a global-scaling normalization method “LogNormalize” 
# that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# In Seurat v5, Normalized values are stored in pbmc[["RNA"]]$data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# if you want these settings, do this with less code:
pbmc <- NormalizeData(pbmc)
# option for different normalization (described in tutorial)

# Identification of highly variable features (feature selection)
# calculate feature that vary wildly between cells in dataset
# helps downstream analysis to discover biological signaling
# directly model the mean-variance relationship inherent in single-cell data

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # finding variable feeatures

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# Scaling the data
# apply a linear transformation ("scaling") to the data
# prior to dimensional reduction techniques like PCA
  # techniques to reduce the dimensions of data through different means (guidelines)
# the ScaleData() function:
  # Shifts the expression of each gene, so that the mean expression across cells is 0
  # Calcs the expression of each gene, so that the variance across cells is 1
    # A standard normal distribution has a mean of 0 and variance of 1. 
    # This is also known as a z distribution.
    # this step gives equal weight in downstream analysis so highly expressed genes don't dominate
  # The results of this are stored in pbmc[["RNA"]]$scale.data
  # By default, only variable features are scaled.
  # You can specify the features argument to scale additional features
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# you have now scaled your data so 


# Perform linear dimensional reduction
# PCA on the scaled data
# default is the highly vairbale features we screened for are used, but we can manually set features if wanted
# For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, 
# representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset.

# Notes on PCA (statquest):
# same concepts as in umap notes: try plotting 4d+ data
# pca takes 4d+ and makes PCA plot to show similar data clusters together, and what variables (genes) account for clustering
# PCA can tell us how accurate a 2d graph is
# imagine an average measurement for gene 1 and 2 (2d data)
# shift points to center their mean on the origin, fit a linear line to it (PC1), and it has a slope
# see the spread of data on the gene 1 axis vs the gene 2 axis (slope tells us this, its like a ratio of gene 1 to gene 2)
# gene 1 is more important
# linear combination = drink ratio of how much the spread of gene 1 vs gene 2 affects their distance from the mean (origin)
# scale the ratio to make PC1 = 1 distance, and PC2 is now scaled down
# unit vector (scaled down ratio) = eigenvector, or singular vector, for PC1
# eigenvalue = sum of squared distances for the best fit line = eigenvalue of PC1
  # just measures of variation
  # proprtion of variation that each PC is accountable for
# e.g. pc1 accoiunts for 80% of total variation around the PCs, and pc2 20%
# PC2 is the line through the origin that is perpendicular to PC1: opposite ratio
# just rotate these two lines to be perpendicular to fine the PCA plot
# process is PCA by SVD
# scree plot is a graphical representation of the percentages of variation that each pc accounts for
  # if skree plot has even spread across many PCs... not good to use 2d PCA
# numbers of PCs is number of variables or number of samples, whichever is smaller
# number of PCs that account for the vast majority of variation would be a good number to use


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# this shows the variance in each treatment group attributable to the PCs
# if the groups have converged, this means that the axis they have converged on is not attributable to their difference
# in this dataset, PC2 could be interpreted as helping the difference between untreated pbmcs, but not for treated pbmcs which have converged
DimPlot(pbmc, reduction = "pca") + NoLegend()

# in particular, DimHeatMap() shows PCs with most heterogenity
# help choose PCs for downstream analysis
# Both cells and features are ordered according to their PCA scores
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, 
# which dramatically speeds plotting for large datasets.
# PC = a parameter calculated on the basis of the variance (pc1 explains most of variance across points, pc2 is second best explaining variance)
# statquest videos can help explains*****
# weight of difference variant factors is different
# one will separate pop into 2, 3, 4 more than another
# percentages - pc1 and 2 account for most variance
# next step could be pathway analysis
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# lets see many PC heat maps now. where does it get messy?
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the ‘dimensionality’ of the dataset
# remove technical noise
# cluster cells based on PCA scores
# combines information across a correlated feature set
# Elbow plot helps us see a ranking of principle components based on the percentage of variance explained by each one
# ‘elbow’ around PC9-10, majority of true signal is captured in the first 10 PCs
ElbowPlot(pbmc)


# Cluster the cells
# distance metric drives clustering analysis based on previously identified PCs
# cluster cells into "Quasi-clusters" or "communities"

# KNN graph = a graph in which two vertices p and q are connected by an edge, 
  # if the distance between p and q is among the k-th smallest distances from p to other objects from P
# construct a KNN graph based on the euclidean distance in PCA space
# FindNeighbors() function, takes as input the previously defined dimensionality of the dataset (first 10 PCs)
# cluster cells together with algorithm, increasing the number of clusters with certain parameters
# resolution parameter =  0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 
# The clusters can be found using the Idents() function.
pbmc <- FindNeighbors(pbmc, dims = 1:10) # low neighbors = small, separated clusters; high neighbors = retins global structure
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# Run non-linear dimensional reduction (UMAP/tSNE)
# visualize and explore these datasets.
# learn underlying structure in the dataset
# place similar cells together in low-dimensional space
# graph-based clusters determined above should co-localize on these dimension reduction plots

# limitations: preserves local distances in the dataset 
# (i.e. ensuring that cells with very similar gene expression profiles co-localize), 
# but often do not preserve more global relationships
# shows cells with closer relationships closer on a 2d or 3d plot
# it moves the projected points around step by step 
  # until the similarities between points in the low-dimensional dataset !!
  # resemble the similarities between cells in the original dataset!!

# UMAP notes (statquest):
# imagine tyring to draw a 4+ dimensional graph... not possible
# PCA only good with first two PCs account for most of the variability of the data
# enter UMAP; take high dimensional data and outputs a low dimensional graph
# umap calcs similarity scores to identify clustered points and preserve is in low dimensional graphs
# umap calcs distance between each pair of high dimensional points; all distances, imagine this on a 1d graph, 
  # and a curve over this representing similarity scores
# e.g. 10 dimensional neighbors, and this includes the point itself
# how does the number of nearest neighbors affect things?
# log2(num neighbors) = curve over 1d plot = defines the curve shape
# y-axis coordinat4es for nearest neighbors add up to the log2(num neighbors)
# repeat for all other points
# b to c has different similarity scofe than c to b, bc sim of similarity scores is always equal to log2(num neighbors)
  # and scaling curves ensures every point is similar to at lease one other point in the dataset
# therefore - > low dimensional graph becomes high dimensional
# umap pciks points to be closer together based on the weights of their similarity scores
# and picks points that should be farther from other points
# how much to move a point? umap calcs low dimensional similarity scores, which come from t-distribution curve (similar bell shape)
# minimize low dimensionality score for far away points, and maximize for close points
# note: moving a closer to b does make it closer to d (another cluster point), but not much relatively

# to compare, t-SNE starts with a random initialization of a low-dimensional graph
# umap uses spectral embedding to initialize the low dimensional graph, every time
# umap better for larger datasets

# umap low number of enighbors is like seeing details, but not big picture
# high num neighbors gives you more of the big picture and less details (more global relationships)

pbmc <- FindNeighbors(pbmc, dims = 1:10) # low neighbors = small, separated clusters; high neighbors = retins global structure
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# You can save the object at this point so that it can easily be loaded back in 
# without having to rerun the computationally intensive steps performed above, 
# or easily shared with collaborators.
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds") #?


# Finding differentially expressed features (cluster biomarkers)
# find markers that define clusters via differential expression (DE)
# positive and negative markers of a single cluster
#  Seurat v5, we use the presto package top improve speed for large datasets

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#ROC test returns the ‘classification power’ for any individual marker 
# (ranging from 0 - random, to 1 - perfect)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0,
                                logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)

# VlnPlot() shows expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()












