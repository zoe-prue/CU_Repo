# ZP 
# PCA on corrected logCPM
# 3/4/25

#### redoing code - haley's code at bottom
#### PCA 1 and 2
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggfortify)

setwd("~/Desktop/CU_coding/RNA-seq/RNA-seq_2/meta_data/")
results_dir <- "~/Desktop/CU_coding/RNA-seq_2/PCA_plots/"

# Read in metadata
#change working direct manually if needed
meta_data <- read.csv("RNAseq_study_design.csv", header = TRUE)
rownames(meta_data) <- meta_data$sample

# Remove sample from meta_data
meta_data_filtered <- meta_data %>%
  filter(sample != "Ep26-JAK1-WT-3")

# Convert the "Sample" column to row names using base R
rownames(meta_data_filtered) <- meta_data_filtered$sample
meta_data_filtered$sample <- NULL  # Remove the redundant column

setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/log2cpm_expression")

# Read in expression data
# can just go to file and setwd
exp_data <- read.table("GE_matrix_corrected_experiment_no_batch_effect_logCPM_voomed.txt", header = TRUE, sep = ",")
colnames(exp_data) <- gsub("\\.", "-", colnames(exp_data))

# Filter exp_data to match meta_data
exp_data_filtered <- exp_data[, colnames(exp_data) %in% rownames(meta_data_filtered)]

# Ensure 'group' and 'experiment' are factors
meta_data_filtered$group <- factor(meta_data_filtered$group)
meta_data_filtered$experiment <- factor(meta_data_filtered$experiment)

# Verify levels
print(levels(meta_data_filtered$group))
print(levels(meta_data_filtered$experiment))

# Check column names match row names
if (!all(colnames(exp_data_filtered) %in% rownames(meta_data_filtered))) {
  stop("Column names in expression data do not match row names in metadata.")
}

# Perform PCA
pca <- prcomp(t(exp_data_filtered), scale. = FALSE)

##### set your wd to PCA_plots manually

# Create PCA plots
pdf(paste0(results_dir, "PCA_no_batch_1_2_scale_FALSE.pdf"), width = 7, height = 5)  # Portrait orientation

# Create PCA plot
p_combined_corrected <- autoplot(pca, data = meta_data_filtered,
                                 colour = 'group', 
                                 shape = 'experiment',
                                 size = 2,
                                 alpha = 0.6, 
                                 x = 1,
                                 y = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right", legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) +
  scale_fill_manual(values = c("#aba1f9", "blue", "#228B22")) +
  ggtitle("PCA 1 and 2 by Group and Experiment no Batch Effect (Scale = False)")

# Print the plot
print(p_combined_corrected)

# Close the PDF
dev.off()

# If you'd prefer to save as PNG, you can do it similarly
png(paste0(results_dir, "PCA_combined_corrected_plot.png"), width = 700, height = 900)  # Adjust size in pixels
print(p_combined)
dev.off()

#### PCA 2 and 3
# Ensure required libraries are loaded
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggfortify)

# Set working directory and results directory
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/")
results_dir <- "~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/PCA/"

# Read in metadata
meta_data <- read.csv("RNAseq_study_design.csv", header = TRUE)
rownames(meta_data) <- meta_data$sample

# Remove sample from meta_data
meta_data_filtered <- meta_data %>%
  filter(sample != "Ep26-JAK1-WT-3")

# Convert the "Sample" column to row names using base R
rownames(meta_data_filtered) <- meta_data_filtered$sample
meta_data_filtered$sample <- NULL  # Remove the redundant column

setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/log2cpm_expression")

# Read in expression data
exp_data <- read.table("GE_matrix_corrected_Group_logCPM_voomed.txt", header = TRUE, sep = ",")
colnames(exp_data) <- gsub("\\.", "-", colnames(exp_data))

# Filter exp_data to match meta_data
exp_data_filtered <- exp_data[, colnames(exp_data) %in% rownames(meta_data_filtered)]

# Ensure 'group' and 'experiment' are factors
meta_data_filtered$group <- factor(meta_data_filtered$group)
meta_data_filtered$experiment <- factor(meta_data_filtered$experiment)

# Verify levels
print(levels(meta_data_filtered$group))
print(levels(meta_data_filtered$experiment))

# Check column names match row names
if (!all(colnames(exp_data_filtered) %in% rownames(meta_data_filtered))) {
  stop("Column names in expression data do not match row names in metadata.")
}

# Perform PCA (transposed to match expression data correctly)
pca <- prcomp(t(exp_data_filtered), scale. = TRUE)

# Create PCA plots for PC2 and PC3
pdf(paste0(results_dir, "PCA_PC2_PC3_corrected.pdf"), width = 7, height = 5)  # Portrait orientation

# Create PCA plot for PC2 and PC3
p_combined_pc2_pc3 <- autoplot(pca, data = meta_data_filtered,
                               colour = 'group', 
                               shape = 'experiment', 
                               size = 2, 
                               alpha = 0.6, 
                               x = 2, # Second principal component (PC2) on the x-axis
                               y = 3) + # Third principal component (PC3) on the y-axis
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right", legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) +
  scale_fill_manual(values = c("#aba1f9", "blue", "#228B22")) +
  ggtitle("PCA by Group and Experiment (PC2 vs PC3)")

# Print the plot
print(p_combined_pc2_pc3)

# Close the PDF
dev.off()

# Save the plot as PNG as well
png(paste0(results_dir, "PCA_PC2_PC3_corrected_plot.png"), width = 700, height = 900)  # Adjust size in pixels
print(p_combined_pc2_pc3)
dev.off()





