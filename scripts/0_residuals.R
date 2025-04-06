library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(utils)
library(edgeR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(limma)

setwd("~/Desktop/CU_coding/RNA-seq/")

# uncorrected counts table w/out outlier
# removing batch effect to get residuals (correcting for batch)
GE <- read.table("DGE_RNAseq/counts_kallisto/Homo_sapiens.GRCh38.87.kallisto.gene_level.lengthScaledTPM_counts_HGNCsymbols.txt", header = TRUE, sep = ",", check.names = FALSE)
colnames(GE) <- gsub("\\.","-",colnames(GE))
meta_data <- read.csv("meta_data/RNAseq_study_design.csv", header = TRUE) # my design
# also factor experiment
meta_data$group = factor(meta_data$group, levels = c("WT","severe_GoF", "mild_GoF"))
meta_data$experiment = factor(meta_data$experiment, levels = c("21","25", "26"))

# ## add scale M_aligned
# meta_data$M_aligned_scale <- scale(meta_data$M_aligned)

# ## subset on protein coding genes
# coding_ids = read.table("genelists/Mouse.mm10.81.HUGOgene_transcriptAssoc.protein_coding.HUGO.txt", header = TRUE)
# GE <- GE[which(rownames(GE) %in% coding_ids$gene_name), ]

## reorder colnames
GE <- GE[meta_data$sample]
colnames(GE) <- meta_data$sample

### remove lowly expressed genes (joaquin method), median expression > 1
dge <- DGEList(counts = GE)
dge <- calcNormFactors(dge)
design <- model.matrix(~ experiment, data = meta_data)
v <- voom(dge, design, plot = TRUE)

tab = data.frame(genes = rownames(GE), medians=apply(v$E, 1, median), order = 1:nrow(GE))
tab = tab[order(-tab$medians), ]
tab$order_by_median = 1:nrow(tab)

## threshold at median = 1
tab = tab[order(tab$order), ]

length(which(rownames(GE) != rownames(tab)))
GE <- GE[which(tab$medians > 1), ]
## 12451 genes

## if you want to save voomed reads with batch effect still intact
# write.table(corrected_expression, "DGE_RNAseq/log2cpm_expression/GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt", quote = FALSE, sep = ",")

## voom after removal of lowly expressed genes
dge <- DGEList(counts = GE)
dge <- calcNormFactors(dge)

## correct for M_aligned_scale
design <- model.matrix(~ experiment, data = meta_data) # design matrix; replace m_aligned_scale whatever experiment batch column
v <- voom(dge, design, plot = TRUE) # fitting a model
vfit <-lmFit(v, design)
vfit <- eBayes(vfit) # beta coefficientss are output of linear model
corrected_expression <- v$E - vfit$coefficients[,"experiment25"]%*%t(design[,"experiment25"]) -vfit$coefficients[,"experiment26"]%*%t(design[,"experiment26"])
write.table(corrected_expression, "DGE_RNAseq/log2cpm_expression/GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt", quote = FALSE, sep = ",")

weights <- v$weights
colnames(weights) <- colnames(corrected_expression)
rownames(weights) <- rownames(corrected_expression)
write.table(weights, paste0("DGE_RNAseq/log2cpm_expression/GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt"), quote = FALSE, sep = ",")

