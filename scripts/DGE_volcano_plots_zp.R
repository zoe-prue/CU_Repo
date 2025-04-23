# Volcano plots 
# DGE RNA-seq 
# ZP 
# 4/22/23

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

devtools::install_github('kevinblighe/EnhancedVolcano')

library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(ggrepel)
library(EnhancedVolcano)

##
setwd("~/Desktop/CU_coding/RNA-seq_2/DEgenes_plots")
results_dir <- "~/Desktop/CU_coding/RNA-seq_2/DEgenes"

DEgenes <- read.csv("DEgenes_by_experiment_comparisons_fdr_0.2.csv")

# Assuming your data is in a dataframe called 'de_data'
volcano_plot <- ggplot(DEgenes, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  facet_wrap(~ comparison, scales = "free") +
  labs(x = "Log2 Fold Change", 
       y = "-Log10(p-value)",
       title = "Facetted Volcano Plots of Differentially Expressed Genes",
       color = "Significant") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, face = "bold"))

print(volcano_plot)

ggsave("volcano_plot.png", 
       plot = volcano_plot,
       width = 12, height = 5, dpi = 300)

# Identify top significant genes for each comparison
top_genes <- DEgenes %>%
  group_by(comparison) %>%
  filter(significant == TRUE) %>%
  arrange(pvalue) %>%
  slice_head(n = 5)  # Top 5 most significant per comparison

# Create the plot with labels
volcano_plot_labeled <- ggplot(DEgenes, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text(data = top_genes, 
            aes(label = hgnc_symbol), 
            size = 3, vjust = 1.5, hjust = 0.5, color = "black") +
  facet_wrap(~ comparison, scales = "free") +
  labs(x = "Log2 Fold Change", 
       y = "-Log10(p-value)",
       title = "Volcano Plots by Comparison with Batch Effect",
       color = "Significant") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, face = "bold"))

print(volcano_plot_labeled)

ggsave("volcano_plot_labeled.png", 
       plot = volcano_plot_labeled,
       width = 12, height = 5, dpi = 300)



