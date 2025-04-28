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

setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq")
results_dir <- "lowly_expressed_genes_removed/"
GE <- read.table("~/Desktop/RNA-seq_big_files/counts_kallistoHomo_sapiens.GRCh38.87.kallisto.gene_level.lengthScaledTPM_counts_HGNCsymbols.txt", header = TRUE, sep = ",")
setwd("~/Desktop/CU_coding/RNA-seq/meta_data")
meta_data <- read.csv("RNAseq_study_design.csv") %>%
  filter(sample != "Ep26-JAK1-WT-3")
meta_data$sample <- gsub("-", ".", meta_data$sample)  # Replace "-" with "."
meta_data$group = factor(meta_data$group, levels = c("WT","mild_GoF", "severe_GoF"))
GE <- GE[, colnames(GE) != "Ep26.JAK1.WT.3"]

## subset on protein coding genes
# coding_ids = read.table("genelists/Mouse.mm10.81.HUGOgene_transcriptAssoc.protein_coding.HUGO.txt", header = TRUE)
# GE <- GE[which(rownames(GE) %in% hgnc_symbol$gene_name), ]

# reorder columns
# check that labels are the same
# colnames(GE)
# meta_data$sample

# GE <- GE[, meta_data$sample, drop = FALSE]  # Reorder columns
# # GE <- GE[meta_data$sample]
# colnames(GE) <- meta_data$sample

# Fix: separate gene names and keep only numeric counts
gene_names <- GE$hgnc_symbol
GE_counts <- GE[, -which(colnames(GE) == "hgnc_symbol")]

### remove lowly expressed genes (joaquin method), median expression > 1
dge <- DGEList(counts = GE_counts)
dge$genes <- data.frame(hgnc_symbol = gene_names)  # reattach gene names
dge <- calcNormFactors(dge)
design <- model.matrix(~ group, data = meta_data)

setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq")
pdf(paste0(results_dir,"voom_before_lowly_expressed_removal.pdf"))
v <- voom(dge, design, plot = TRUE)
dev.off()

tab = data.frame(genes = GE$hgnc_symbol, medians=apply(v$E, 1, median), order = 1:nrow(GE))
tab = tab[order(-tab$medians), ]
tab$order_by_median = 1:nrow(tab)

pdf(paste0(results_dir,"lowly_expressed_genes.pdf"))
plot(tab$order_by_median, tab$medians)
dev.off()
## threshold at median = 1
tab = tab[order(tab$order), ]

length(which(rownames(GE) != rownames(tab)))
GE <- GE[which(tab$medians > 1), ]
## 11959 genes

# Fix: separate gene names and keep only numeric counts
gene_names <- GE$hgnc_symbol
GE_counts <- GE[, -which(colnames(GE) == "hgnc_symbol")]

### remove lowly expressed genes (joaquin method), median expression > 1
dge <- DGEList(counts = GE_counts)
dge$genes <- data.frame(hgnc_symbol = gene_names)  # reattach gene names
dge <- calcNormFactors(dge)
design <- model.matrix(~ group, data = meta_data)

pdf(paste0(results_dir,"voom_after_lowly_expressed_removal.pdf"))
v <- voom(dge, design, plot = TRUE)
dev.off()
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq")



write.table(v, "log2cpm_expression/GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt", quote = FALSE, sep = ",")


### STOP HERE FOR NOW, PLOT PCA ON UNCORRECTED VALUES ###



# Ensure the design matrix and coefficients are properly matched
design_simple <- model.matrix(~ group, data = meta_data)  # Create the design matrix for group

# Perform voom transformation on the data
v <- voom(dge, design_simple, plot = TRUE)

# Fit the linear model and apply eBayes
vfit <- lmFit(v, design_simple)
vfit <- eBayes(vfit)

# Correct for the group effects (using the relevant coefficients)
# Let's assume you want to remove the effect of the 'group' variable, which is encoded by the columns: 'groupmild_GoF' and 'groupsevere_GoF'
# 1. Subtract group effects
corrected_expression <- v$E - vfit$coefficients[,2] %*% t(design_simple[,2]) -
  vfit$coefficients[,3] %*% t(design_simple[,3])

# 2. Add hgnc_symbol column
corrected_df <- cbind(hgnc_symbol = v$genes$hgnc_symbol, corrected_expression)

# 3. Save to file
write.table(corrected_df, file = "GE_matrix_corrected_experiment_no_batch_effect_logCPM_voomed_2.txt",
            quote = FALSE, sep = ",", row.names = FALSE)

