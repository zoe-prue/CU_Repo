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

# Combine batch and nobatch data frames
volcano_combined <- bind_rows(batch_df, nobatch_expanded)

# Combine with batch data
volcano_combined_sig <- bind_rows(batch_df, nobatch_expanded) %>%
  mutate(significant = (abs(logFC) > beta_thresh & pvalue < p_thresh & fdr < fdr_thresh)) %>%
  filter(significant == TRUE)

# Now proceed with your filtering and plotting using volcano_combined
# 
# # Step 1: Filter for significant DEGs
# filtered <- volcano_combined %>%
#   filter(pvalue < 0.05, abs(logFC) > 0.5)

gene_exp_counts <- volcano_combined_sig %>%
  filter(significant == TRUE) %>%       # Keep only significant rows
  distinct(hgnc_symbol, experiment) %>% # Get unique gene-experiment pairs
  group_by(hgnc_symbol) %>%              # Group by gene
  summarise(n_exp = n())

# Step 3: Get genes significant in >= 2 experiments
genes_in_2_exp <- gene_exp_counts %>%
  filter(n_exp >= 2) %>%
  pull(hgnc_symbol)

# Step 4: Filter original data to include only those genes
filtered_2exp <- volcano_combined_sig %>%
  filter(hgnc_symbol %in% genes_in_2_exp, source == "Individual Experiment")

# Make sure significant is logical
# volcano_combined <- volcano_combined %>% mutate(significant = as.logical(significant))
# 
# # Count number of distinct experiments per gene where significant == TRUE
# sig_exp_counts <- filtered %>%
#   filter(significant) %>%
#   distinct(hgnc_symbol, experiment) %>%
#   group_by(hgnc_symbol) %>%
#   summarise(n_exp = n())
# 

# Combine with batch data
# volcano_combined <- bind_rows(batch_df, nobatch_expanded) %>%
  # mutate(significant = (abs(logFC) > beta_thresh & pvalue < p_thresh))

# # Step 1: Get genes significant in exactly 2 experiments
# genes_in_2_exp <- volcano_combined %>%
#   filter(significant) %>%
#   distinct(hgnc_symbol, experiment) %>%
#   group_by(hgnc_symbol) %>%
#   summarise(n_exp = n()) %>%
#   filter(n_exp == 2) %>%
#   pull(hgnc_symbol)

library(stringr)

volcano_combined_label <- filtered_2exp %>%
  filter(significant, hgnc_symbol %in% genes_in_2_exp) %>%
  mutate(
    hgnc_symbol = str_trim(hgnc_symbol),
    experiment = str_trim(experiment),
    label = paste0(hgnc_symbol, " (", experiment, ")")
  ) %>%
  group_by(hgnc_symbol) %>%
  filter(n() > 1) %>%
  ungroup()


# Step 3: Plot
p1 <- ggplot(volcano_combined_label, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(color = "red", alpha = 0.7) +
  geom_text_repel(aes(label = label),
                  size = 5,
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Significant Genes in >= 2 Experiments",
       x = "log Fold Change (logFC)",
       y = "-log10(p-value)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 26, hjust = 0.5)# Make facet titles bigger and bold
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p1)

p1 <- p1 + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "white", color = NA)
)


ggsave(
  "volcano_common_DEGenes_dots_all_exps.png",
  plot = p1,
  width = 14,
  height = 10,
  dpi = 300  # Optional: set resolution for high-quality output
)


# # Ensure labels are clean and consistent
# volcano_combined_comparisons <- filtered_2exp %>%
#   filter(significant, hgnc_symbol %in% genes_in_2_exp) %>%
#   mutate(
#     hgnc_symbol = str_trim(hgnc_symbol),
#     experiment = str_trim(experiment),
#     comparison = str_trim(comparison),  # ensure comparison is clean
#     label = paste0(hgnc_symbol, " (", experiment, ")")
#   )

# Plot with faceting by comparison
p2 <- ggplot(volcano_combined_label, aes(x = logFC, y = -log10(pvalue))) +
  geom_point(color = "red", alpha = 0.7) +
  facet_wrap(~ comparison, scales = "free_x", ncol=2) +
  geom_text_repel(
    aes(label = label),
    size = 5,
    na.rm = TRUE,
    max.overlaps = Inf
  ) +
  facet_wrap(~ comparison, scales = "free") +  # Facet by comparison
  theme_minimal() +
  labs(
    title = "Significant Genes in >= 2 Experiments, Facetted by Comparison",
    x = "log Fold Change (logFC)",
    y = "-log10(p-value)"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 26, hjust = 0.5)# Make facet titles bigger and bold
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(p2)

p2 <- p2 + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  strip.background = element_rect(fill = "white", color = NA)
)


ggsave(
  "volcano_common_DEGenes_dots_all_exps_w_comparisons.png",
  plot = p2,
  width = 14,
  height = 10,
  dpi = 300  # Optional: set resolution for high-quality output
)


