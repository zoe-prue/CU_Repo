# Simplified Differential Gene Expression Analysis
library(limma)
library(tidyverse)

setwd("/Users/zoeprue/Desktop/CU_coding/RNA-seq_2/")

# Set parameters
fdr_threshold <- 0.2
results_dir <- "DEgenes/"

# Load data - manually have to setwd for each one
residuals <- read.table("log2cpm_expression/GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt", 
                       header = TRUE, sep = ",", row.names = 1)
rownames(residuals) <- residuals$hgnc_symbol
residuals$hgnc_symbol <- NULL

summary(rowMeans(residuals))

meta_data <- read.csv("meta_data/RNAseq_study_design.csv") %>%
  mutate(group = factor(group, levels = c("WT", "mild_GoF", "severe_GoF"))) %>%
  column_to_rownames("sample")

# Clean sample names
clean_names <- function(x) gsub("^Ep(\\d+)", "Exp\\1", gsub("\\.", "-", x))
colnames(residuals) <- clean_names(colnames(residuals))
rownames(meta_data) <- clean_names(rownames(meta_data))

# Remove problematic sample
residuals <- residuals[, !colnames(residuals) %in% "Ep26-JAK1-WT-3"]
meta_data <- meta_data[!rownames(meta_data) %in% "Exp26-JAK1-WT-3",]

# Define comparisons
comparisons <- list(
  WT_vs_severe = c("WT", "severe_GoF"),
  WT_vs_mild = c("WT", "mild_GoF"),
  mild_vs_severe = c("mild_GoF", "severe_GoF")
)

# Main analysis function
run_dge <- function(comp_name, groups) {
  # Subset metadata
  meta_comp <- meta_data %>% 
    filter(group %in% groups) %>%
    mutate(group = factor(group, levels = groups))
  
  # Skip if not enough samples
  if (n_distinct(meta_comp$group) < 2) {
    message("Skipping ", comp_name, " - insufficient groups")
    return(NULL)
  }
  
  # Align samples
  samples <- intersect(rownames(meta_comp), colnames(residuals))
  expr_data <- residuals[, samples]
  meta_comp <- meta_comp[samples, ]
  meta_comp$experiment <- factor(meta_comp$experiment)
  
  # Create design matrix (will auto-name coefficients as groupLEVEL)
  design <- model.matrix(~ experiment + group, data = meta_comp)
  cat("\nDesign matrix for", comp_name, ":\n")
  print(head(design))
  
  # Fit model - coefficient will be "groupsevere_GoF" etc.
  fit <- lmFit(expr_data, design) %>% eBayes()
  coef_name <- paste0("group",groups[2])  # Gets "groupsevere_GoF" format
  
  # Get results
  results <- tibble(
    gene_id = rownames(fit$coefficients),
    logFC = fit$coefficients[, coef_name],
    pvalue = fit$p.value[, coef_name],
    fdr = p.adjust(fit$p.value[, coef_name], method = "BH"),
    significant = fdr < fdr_threshold,
    comparison = comp_name
  )
}

# Run all comparisons
results <- map2_dfr(names(comparisons), comparisons, run_dge)

# Save results
write_csv(results, file.path(results_dir, paste0("DEgenes_no_batch_effect_all_comparisons_", fdr_threshold, ".csv")))

# Count significant genes
sig_counts <- results %>%
  group_by(comparison) %>%
  summarize(significant_genes = sum(significant, na.rm = TRUE))

print(sig_counts)


