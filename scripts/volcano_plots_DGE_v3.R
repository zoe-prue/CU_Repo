# volcano plots DGE zp v3
# 4/30.25
# ZP

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
p_thresh = 0.05

# setwd
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/DEgenes/")
degenes_in_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenesDGE_results_all_comparisons_fdr_0.2.csv")
degenes_in_batches$gene_id <- NULL
head(degenes_in_batches)
degenes_no_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenes_no_batch_effect_all_comparisons_w_nonsig0.2.csv")
head(degenes_no_batches)
# make column names match
names(degenes_no_batches)[names(degenes_no_batches) == "gene_id"] <- "hgnc_symbol"

# Prep batch data (includes experiments)
batch_df <- degenes_in_batches %>%
  select(hgnc_symbol, comparison, experiment, logFC, fdr, pvalue) %>%
  mutate(source = "batch")  # tag as batch

# Prep no-batch data (same values repeated across experiments)
nobatch_df <- degenes_no_batches %>%
  select(hgnc_symbol, comparison, logFC, fdr, pvalue) %>%
  mutate(source = "no_batch")

# Combine batch and no-batch
volcano_combined <- bind_rows(batch_df, nobatch_df)

# Add significance flag
volcano_combined_sig <- volcano_combined %>%
  mutate(significant = (abs(logFC) > beta_thresh & pvalue < p_thresh))

# Get unique experiment values from batch data
experiments <- unique(batch_df$experiment)

# Expand no_batch data to all experiments
nobatch_expanded <- nobatch_df %>%
  tidyr::crossing(experiment = experiments)  # replicate across experiments

# Combine with batch data
volcano_combined <- bind_rows(batch_df, nobatch_expanded) %>%
  mutate(significant = (abs(logFC) > beta_thresh & fdr < fdr_thresh))

head(volcano_combined)

ggplot(volcano_combined, aes(x = logFC, y = -log10(pvalue),
                             color = source,
                             alpha = significant)) +
  geom_point() +
  facet_wrap(~ experiment) +
  geom_vline(xintercept = -log10(p_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  geom_text_repel(data = filter(volcano_combined, significant),
                  aes(label = hgnc_symbol),
                  size = 2.5,
                  max.overlaps = 100) +
  scale_color_manual(values = c("batch" = "darkred", "no_batch" = "steelblue")) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  labs(title = "Volcano Plot: Batch vs No Batch by Experiment",
       x = "log2 Fold Change", y = "-log10 p-value") +
  theme_minimal()


# Save
ggsave("pvalue0.05_logfc0.5_batch_vs_nobatch_facet.pdf", width = 12, height = 6)


