# Differential Gene Expression Analysis Plots
# ZP
# Updated: May 2024

#######################
## LIBRARY IMPORTS ###
#######################
library(tidyverse)
library(ggrepel)
library(readr)

#######################
## DATA PREPARATION ###
#######################

# Set global parameters
FDR_THRESHOLD <- 0.2
BETA_THRESHOLD <- 0.5
BASE_DIR <- "~/Desktop/CU_coding/RNA-seq_2"

# Read metadata
meta_data <- read_csv(file.path(BASE_DIR, "meta_data/RNAseq_study_design.csv")) %>%
  mutate(group = factor(group, levels = c("WT", "mild_GoF", "severe_GoF"))) %>%
  column_to_rownames("sample")

# Read DE genes data
degenes_batch <- read_csv(file.path(BASE_DIR, "DEgenes/DEgenes_by_experiment_comparisons_fdr_0.2.csv"))
degenes_nobatch <- read_csv(file.path(BASE_DIR, "DEgenes/DEgenes_no_batch_effect_all_comparisons_0.2.csv")) %>%
  rename(hgnc_symbol = gene_id)

# Read and process non-DE genes data
read_expression_matrix <- function(path) {
  read_csv(path) %>%
    select(hgnc_symbol, everything()) %>%
    pivot_longer(
      cols = -hgnc_symbol,
      names_to = "sample",
      values_to = "expression"
    ) %>%
    left_join(
      meta_data %>% rownames_to_column("sample") %>% select(sample, group),
      by = "sample"
    )
}

de_batch <- degenes_batch %>%
  filter(comparison == comparison_name) %>%
  select(hgnc_symbol, pvalue, logFC, experiment) %>%
  mutate(gene_type = "DE")

nonde_batch_sub <- nonde_batch %>%
  filter(comparison == comparison_name) %>%
  filter(!hgnc_symbol %in% de_batch$hgnc_symbol) %>%
  select(hgnc_symbol, pvalue, experiment) %>%
  mutate(gene_type = "non-DE")

de_nobatch <- degenes_nobatch %>%
  filter(comparison == comparison_name) %>%
  select(hgnc_symbol, pvalue, logFC) %>%
  mutate(gene_type = "DE")

nonde_nobatch_sub <- nonde_nobatch %>%
  filter(comparison == comparison_name) %>%
  filter(!hgnc_symbol %in% de_nobatch$hgnc_symbol) %>%
  select(hgnc_symbol, pvalue) %>%
  mutate(gene_type = "non-DE")



#############################
## DATA PROCESSING ###
#############################

# Function to prepare comparison data
prepare_comparison_data <- function(comparison_name) {
  
  # comparison_name <- "WT_vs_severe"
  # DE genes
  de_batch <- degenes_batch %>%
    filter(comparison == comparison_name) %>%
    select(hgnc_symbol, pvalue, logFC, experiment) %>%
    mutate(gene_type = "DE")
  
  de_nobatch <- degenes_nobatch %>%
    filter(comparison == comparison_name) %>%
    select(hgnc_symbol, pvalue, logFC) %>%
    mutate(gene_type = "DE")
  
  # Non-DE genes (assuming all genes not in DE lists are non-DE)
  # This would need adjustment based on your actual data structure
  nonde_batch_sub <- nonde_batch %>%
    filter(!hgnc_symbol %in% de_batch$hgnc_symbol) %>%
    select(hgnc_symbol, pvalue, experiment) %>%
    mutate(gene_type = "non-DE")
  
  nonde_nobatch_sub <- nonde_nobatch %>%
    filter(!hgnc_symbol %in% de_nobatch$hgnc_symbol) %>%
    select(hgnc_symbol, pvalue) %>%
    mutate(gene_type = "non-DE")
  
  # Merge batch and no-batch data
  merged_data <- full_join(
    bind_rows(de_batch, nonde_batch_sub),
    bind_rows(de_nobatch, nonde_nobatch_sub) %>% rename(pvalue_nobatch = pvalue),
    by = "hgnc_symbol"
  )
  
  return(merged_data)
}

# Prepare data for each comparison
wt_severe_data <- prepare_comparison_data("WT_vs_severe")
wt_mild_data <- prepare_comparison_data("WT_vs_mild")
mild_severe_data <- prepare_comparison_data("mild_vs_severe")

#############################
## PLOTTING FUNCTIONS ###
#############################

# Function to create p-value comparison plot
create_pvalue_plot <- function(data, comparison_name) {
  ggplot(data, aes(x = -log10(pvalue), y = -log10(pvalue_nobatch))) +
    # Plot non-DE genes first (as background)
    geom_point(data = filter(data, gene_type == "non-DE"), 
               color = "red", alpha = 0.2, size = 1.5) +
    # Plot DE genes on top
    geom_point(data = filter(data, gene_type == "DE"), 
               color = "black", alpha = 0.7, size = 2) +
    # Label only DE genes
    geom_text_repel(
      data = filter(data, gene_type == "DE"),
      aes(label = hgnc_symbol),
      size = 3,
      color = "black",
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    facet_wrap(~ experiment, ncol = 2) +
    labs(
      x = "-log10 p-value (with batch effect)",
      y = "-log10 p-value (no batch effect)",
      title = paste0(comparison_name, ": Batch vs No Batch Effect"),
      subtitle = "DE genes (black), non-DE genes (red)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# Function to create logFC comparison plot
create_logfc_plot <- function(data, comparison_name) {
  # This would use similar structure but plot logFC instead of p-values
  # Implementation would depend on your specific logFC data
}

#############################
## GENERATE AND SAVE PLOTS ##
#############################

# Create and save p-value plots
p_wt_severe <- create_pvalue_plot(wt_severe_data, "WT vs Severe")
# ggsave("wt_vs_severe_pvalue_comparison.pdf", p_wt_severe, width = 10, height = 8)

p_wt_mild <- create_pvalue_plot(wt_mild_data, "WT vs Mild")
ggsave("wt_vs_mild_pvalue_comparison.pdf", p_wt_mild, width = 10, height = 8)

p_mild_severe <- create_pvalue_plot(mild_severe_data, "Mild vs Severe")
ggsave("mild_vs_severe_pvalue_comparison.pdf", p_mild_severe, width = 10, height = 8)

# Similarly create and save logFC plots
# ...


