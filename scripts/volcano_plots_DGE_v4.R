# volcano plots DGE zp v4
# Modified based on feedback
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

## read in data
results_dir <- "~/Desktop/CU_coding/RNA-seq_2/DEgenes"
meta_data <- read.csv("meta_data/RNAseq_study_design.csv") %>%
  mutate(group = factor(group, levels = c("WT", "mild_GoF", "severe_GoF"))) %>%
  column_to_rownames("sample")
fdr_thresh = 0.2
p_thresh = 0.05
beta_thresh = 0.5  # Added based on your save filename

# setwd
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/DEgenes/")
degenes_in_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenesDGE_results_all_comparisons_fdr_0.2.csv")
degenes_in_batches$gene_id <- NULL
degenes_no_batches <- read_csv("~/Desktop/CU_coding/RNA-seq_2/DEgenes/DEgenes_no_batch_effect_all_comparisons_w_nonsig0.2.csv")

# make column names match
names(degenes_no_batches)[names(degenes_no_batches) == "gene_id"] <- "hgnc_symbol"

# Prep batch data (individual experiments)
batch_df <- degenes_in_batches %>%
  select(hgnc_symbol, comparison, experiment, logFC, fdr, pvalue) %>%
  mutate(source = "Individual Experiment")  # Changed label

# Prep no-batch data (combined analysis)
nobatch_df <- degenes_no_batches %>%
  select(hgnc_symbol, comparison, logFC, fdr, pvalue) %>%
  mutate(source = "Combined Analysis")  # Changed label

# Get unique experiment values from batch data
experiments <- unique(batch_df$experiment)

# Expand no_batch data to all experiments for faceting
nobatch_expanded <- nobatch_df %>%
  tidyr::crossing(experiment = experiments)

# Combine with batch data
volcano_combined <- bind_rows(batch_df, nobatch_expanded) %>%
  mutate(significant = (abs(logFC) > beta_thresh & pvalue < p_thresh & fdr < fdr_thresh))

# Create the volcano plot with larger width ##############################
p <- ggplot(volcano_combined, aes(x = logFC, y = -log10(pvalue),
                                  color = source,
                                  alpha = significant)) +
  geom_point(size = 1) +  # Slightly smaller points
  facet_wrap(~ experiment + comparison, ncol = 3) +  # Better organization
  geom_vline(xintercept = c(-beta_thresh, beta_thresh), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "gray40") +
  geom_text_repel(
    data = filter(volcano_combined, significant),
    aes(label = hgnc_symbol),
    size = 2,
    max.overlaps = 100,  # Reduced from 100
    box.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Set x-axis breaks at increments of 1
  scale_x_continuous(breaks = seq(floor(min(volcano_combined$logFC)), 
                                  ceiling(max(volcano_combined$logFC)), 
                                  by = 1)) +
  scale_color_manual(
    values = c("Individual Experiment" = "darkred", 
               "Combined Analysis" = "steelblue")
  ) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.2)) +
  labs(
    title = "Differential Expression Analysis Comparison",
    subtitle = "Individual Experiments vs. Combined Analysis",
    x = "log2 Fold Change", 
    y = "-log10 p-value",
    color = "Analysis Type",
    alpha = "Significant"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

print(p)

# Save with larger width
ggsave(
  "volcano_comparison_individual_vs_combined.pdf", 
  plot = p, 
  width = 14,  # Increased from 12
  height = 10  # Increased from 6
)


volcano_combined_sig <- volcano_combined %>%
  filter(significant == TRUE)

head(volcano_combined_sig)



# Create the volcano plot with non-significant genes as black dots ##################
p2 <- ggplot(volcano_combined, aes(x = logFC, y = -log10(pvalue))) +
  # Plot all non-significant genes as black dots first (background)
  geom_point(data = filter(volcano_combined, !significant), 
             color = "black", alpha = 0.2, size = 1) +
  # Then plot significant genes colored by source
  geom_point(data = filter(volcano_combined, significant), 
             aes(color = source), alpha = 0.8, size = 1) +
  facet_wrap(~ experiment + comparison, ncol = 3) +
  geom_vline(xintercept = c(-beta_thresh, beta_thresh), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "gray40") +
  geom_text_repel(
    data = filter(volcano_combined, significant),
    aes(label = hgnc_symbol, color = source),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = 'grey50'
  ) +
  # Set x-axis breaks at increments of 1
  scale_x_continuous(breaks = seq(floor(min(volcano_combined$logFC)), 
                                  ceiling(max(volcano_combined$logFC)), 
                                  by = 1)) +
  scale_color_manual(
    values = c("Individual Experiment" = "darkred", 
               "Combined Analysis" = "steelblue"),
    name = "Analysis Type"
  ) +
  labs(
    title = "Differential Expression Analysis Comparison",
    subtitle = "Individual Experiments vs. Combined Analysis with Experiment Correction",
    x = "log2 Fold Change", 
    y = "-log10 p-value"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

# # Save with larger width
# ggsave(
#   "volcano_comparison_individual_vs_combined_black_nonsig_dots.pdf", 
#   plot = p2, 
#   width = 14,
#   height = 10
# )

p2 <- p2 + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "white", color = NA)
)

print(p2)

ggsave(
  "volcano_comparison_individual_vs_combined_black_nonsig_dots.png",
  plot = p2,
  width = 14,
  height = 10,
  dpi = 300  # Optional: set resolution for high-quality output
)

