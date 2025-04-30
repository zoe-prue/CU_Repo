# plots for DGE
# ZP
# 4/30/25

## Library calls
library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(ggrepel)
library(readr)
library(tidyverse)
library(ggrepel)  # install.packages("ggrepel") if needed

## read in data
results_dir <- "~/Desktop/CU_coding/RNA-seq_2/DEgenes"
meta_data <- read.csv("meta_data/RNAseq_study_design.csv") %>%
  mutate(group = factor(group, levels = c("WT", "mild_GoF", "severe_GoF"))) %>%
  column_to_rownames("sample")
fdr_thresh = 0.2
beta_thresh = 0.5

# setwd
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/DEgenes/")
degenes_in_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenes_by_experiment_comparisons_fdr_0.2.csv")
degenes_in_batches$gene_id <- NULL
head(degenes_in_batches)
degenes_no_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenes_no_batch_effect_all_comparisons_0.2.csv")
head(degenes_no_batches)
# make column names match
names(degenes_no_batches)[names(degenes_no_batches) == "gene_id"] <- "hgnc_symbol"

############### PVALUES #######################################################

# Subset for WT_vs_severe ################################################
wt_severe_batch <- subset(degenes_in_batches, comparison == "WT_vs_severe")
library(dplyr)
library(ggplot2)

# Subset no-batch data for WT_vs_severe (assumed same across experiments)
wt_severe_nobatch <- degenes_no_batches %>%
  filter(comparison == "WT_vs_severe") %>%
  select(hgnc_symbol, pvalue)

names(wt_severe_nobatch)[2] <- "pvalue_nobatch"  # base R rename

# Subset batch data for WT_vs_severe (has multiple experiments)
wt_severe_batch <- degenes_in_batches %>%
  filter(comparison == "WT_vs_severe") %>%
  select(hgnc_symbol, pvalue, experiment)

names(wt_severe_batch)[2] <- "pvalue_batch"  # base R rename

# Merge batch and no-batch p-values by gene
merged_all <- merge(wt_severe_batch, wt_severe_nobatch, by = "hgnc_symbol")

library(ggplot2)
library(ggrepel)  # install.packages("ggrepel") if needed

ggplot(merged_all, aes(x = -log10(pvalue_batch), y = -log10(pvalue_nobatch), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +  # avoids overlap
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "-log10 p-value (with batch effect)",
    y = "-log10 p-value (no batch effect)",
    title = "WT vs Severe: Batch vs No Batch (by Experiment)"
  ) +
  coord_cartesian(clip = "off") +  # Prevent text cutoff
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)  # top, right, bottom, left
  )

ggsave("wt_vs_severe_batch_vs_nobatch.pdf", width = 10, height = 6)


# Subset for WT_vs_mild ###############################################
# Subset for WT_vs_mild
wt_mild_batch <- subset(degenes_in_batches, comparison == "WT_vs_mild")

library(dplyr)
library(ggplot2)
library(ggrepel)  # install.packages("ggrepel") if needed

# Subset no-batch data for WT_vs_mild (assumed same across experiments)
wt_mild_nobatch <- degenes_no_batches %>%
  filter(comparison == "WT_vs_mild") %>%
  select(hgnc_symbol, pvalue)

names(wt_mild_nobatch)[2] <- "pvalue_nobatch"  # base R rename

# Subset batch data for WT_vs_mild (has multiple experiments)
wt_mild_batch <- degenes_in_batches %>%
  filter(comparison == "WT_vs_mild") %>%
  select(hgnc_symbol, pvalue, experiment)

names(wt_mild_batch)[2] <- "pvalue_batch"  # base R rename

# Merge batch and no-batch p-values by gene
merged_mild <- merge(wt_mild_batch, wt_mild_nobatch, by = "hgnc_symbol")

# Optional: define significance (e.g., FDR < 0.2, or p < 0.05)
merged_mild$significant <- merged_mild$pvalue_nobatch < 0.05  # adjust threshold as needed

# Plot
ggplot(merged_mild, aes(x = -log10(pvalue_batch), y = -log10(pvalue_nobatch), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "-log10 p-value (with batch effect)",
    y = "-log10 p-value (no batch effect)",
    title = "WT vs Mild: Batch vs No Batch (by Experiment)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)
  )

# Save to PDF
ggsave("wt_vs_mild_batch_vs_nobatch.pdf", width = 10, height = 6)

# Subset mild vs severe ################################################

# Subset for mild_vs_severe
mild_severe_batch <- subset(degenes_in_batches, comparison == "mild_vs_severe")

# Subset no-batch data for mild_vs_severe
mild_severe_nobatch <- degenes_no_batches %>%
  filter(comparison == "mild_vs_severe") %>%
  select(hgnc_symbol, pvalue)

names(mild_severe_nobatch)[2] <- "pvalue_nobatch"  # base R rename

# Subset batch data for mild_vs_severe (has multiple experiments)
mild_severe_batch <- degenes_in_batches %>%
  filter(comparison == "mild_vs_severe") %>%
  select(hgnc_symbol, pvalue, experiment)

names(mild_severe_batch)[2] <- "pvalue_batch"  # base R rename

# Merge batch and no-batch p-values by gene
merged_mild_severe <- merge(mild_severe_batch, mild_severe_nobatch, by = "hgnc_symbol")

# Optional: define significance (e.g., p < 0.05 for no-batch p-values)
merged_mild_severe$significant <- merged_mild_severe$pvalue_nobatch < 0.05

# Plot
ggplot(merged_mild_severe, aes(x = -log10(pvalue_batch), y = -log10(pvalue_nobatch), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "-log10 p-value (with batch effect)",
    y = "-log10 p-value (no batch effect)",
    title = "Mild vs Severe: Batch vs No Batch (by Experiment)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)
  )

# Save to PDF
ggsave("mild_vs_severe_batch_vs_nobatch.pdf", width = 10, height = 6)

############### BETAS #######################################################

# Subset no-batch data for WT_vs_mild
wt_mild_nobatch <- degenes_no_batches %>%
  filter(comparison == "WT_vs_mild") %>%
  select(hgnc_symbol, logFC)

names(wt_mild_nobatch)[2] <- "logFC_nobatch"

# Subset batch data for WT_vs_mild (multiple experiments)
wt_mild_batch <- degenes_in_batches %>%
  filter(comparison == "WT_vs_mild") %>%
  select(hgnc_symbol, logFC, experiment)

names(wt_mild_batch)[2] <- "logFC_batch"

# Merge
merged_logFC <- merge(wt_mild_batch, wt_mild_nobatch, by = "hgnc_symbol")

# Plot
ggplot(merged_logFC, aes(x = logFC_batch, y = logFC_nobatch, color = abs(logFC_batch) > 1 | abs(logFC_nobatch) > 1)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "logFC (with batch effect)",
    y = "logFC (no batch effect)",
    title = "WT vs Mild: logFC (Batch vs No Batch)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave("wt_vs_mild_logFC_batch_vs_nobatch.pdf", width = 10, height = 6)


#########

# Subset no-batch data for WT_vs_severe
wt_severe_nobatch <- degenes_no_batches %>%
  filter(comparison == "WT_vs_severe") %>%
  select(hgnc_symbol, logFC)

names(wt_severe_nobatch)[2] <- "logFC_nobatch"

# Subset batch data for WT_vs_severe (multiple experiments)
wt_severe_batch <- degenes_in_batches %>%
  filter(comparison == "WT_vs_severe") %>%
  select(hgnc_symbol, logFC, experiment)

names(wt_severe_batch)[2] <- "logFC_batch"

# Merge
merged_logFC_severe <- merge(wt_severe_batch, wt_severe_nobatch, by = "hgnc_symbol")

# Plot
ggplot(merged_logFC_severe, aes(x = logFC_batch, y = logFC_nobatch, color = abs(logFC_batch) > 1 | abs(logFC_nobatch) > 1)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "logFC (with batch effect)",
    y = "logFC (no batch effect)",
    title = "WT vs Severe: logFC (Batch vs No Batch)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave("wt_vs_severe_logFC_batch_vs_nobatch.pdf", width = 10, height = 6)

###########

# Subset no-batch data for mild_vs_severe
mild_severe_nobatch <- degenes_no_batches %>%
  filter(comparison == "mild_GoF_vs_severe_GoF") %>%
  select(hgnc_symbol, logFC)

names(mild_severe_nobatch)[2] <- "logFC_nobatch"

# Subset batch data for mild_vs_severe
mild_severe_batch <- degenes_in_batches %>%
  filter(comparison == "mild_GoF_vs_severe_GoF") %>%
  select(hgnc_symbol, logFC, experiment)

names(mild_severe_batch)[2] <- "logFC_batch"

# Merge
merged_logFC_mild_severe <- merge(mild_severe_batch, mild_severe_nobatch, by = "hgnc_symbol")

# Plot
ggplot(merged_logFC_mild_severe, aes(x = logFC_batch, y = logFC_nobatch, color = abs(logFC_batch) > 1 | abs(logFC_nobatch) > 1)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = hgnc_symbol), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  facet_wrap(~ experiment, scales = "free_y") +
  labs(
    x = "logFC (with batch effect)",
    y = "logFC (no batch effect)",
    title = "Mild vs Severe: logFC (Batch vs No Batch)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave("mild_vs_severe_logFC_batch_vs_nobatch.pdf", width = 10, height = 6)


