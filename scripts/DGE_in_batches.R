library(limma)
library(dplyr)
library(tidyverse)
library(tibble)

###### Variables to change ########
fdr_threshold <- 0.2
results_dir <- "~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/DGE_outputs/DGE_by_experiment_comparisons/"
###### End variables #############

# Load and prepare data ----
# Read residuals (log2CPM)
setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/log2cpm_expression/")
residuals <- read.table("GE_matrix_corrected_experiment_with_batch_effect_logCPM_voomed.txt", 
                        header = TRUE, sep = ",", row.names = 1)

# Read metadata
# factor the groups and experiments
# tmake the sample column actually into the rownames of the meta_data -> easier for later :)
setwd("~/Desktop/CU_coding/RNA-seq/meta_data")
meta_data <- read.csv("RNAseq_study_design.csv") %>%
  mutate(
    experiment = factor(experiment, levels = c(21, 25, 26)),
    group = factor(group, levels = c("WT", "mild_GoF", "severe_GoF"))
  ) %>%
  column_to_rownames("sample")


# Clean sample names ----
colnames(residuals) <- gsub("\\.", "-", colnames(residuals))

# Remove problematic sample if exists from residuals
if ("Ep26-JAK1-WT-3" %in% colnames(residuals)) {
  residuals <- residuals[, !colnames(residuals) %in% "Ep26-JAK1-WT-3"]
}
"Ep26-JAK1-WT-3" %in% colnames(residuals)  # Should return FALSE

# Remove the bad sample by name from meta_data
meta_data <- meta_data[rownames(meta_data) != "Ep26-JAK1-WT-3", ]
"Ep26-JAK1-WT-3" %in% rownames(meta_data)  # Should return FALSE

# Convert ALL "EpXX-" to "ExpXX-" in residuals column names
colnames(residuals) <- gsub("^Ep(\\d+)", "Exp\\1", colnames(residuals))
# Verify changes
grep("^Ep", colnames(residuals), value = TRUE)  # Should return character(0)
grep("^Exp", colnames(residuals), value = TRUE)  # Should show Exp26 samples

# Convert ALL "EpXX-" to "ExpXX-" in meta_data column names
rownames(meta_data) <- gsub("^Ep(\\d+)", "Exp\\1", rownames(meta_data))
# Verify changes
grep("^Ep", rownames(meta_data), value = TRUE)  # Should return character(0)
grep("^Exp", rownames(meta_data), value = TRUE)  # Should show Exp26 samples

# Main analysis loop ----
all_results <- list()
significant_counts <- list()

# first subset by experiment
for (exp in meta_data$experiment) {
  cat("\n=== Processing experiment", exp, "=== :-) \n")
  
  # Subset data
  meta_exp <- meta_data %>%
    filter(experiment == exp) %>%
    mutate(group = droplevels(group)) # dropping other unwanted experiments
  
  # Check group representation
  group_counts <- table(meta_exp$group)
  print(group_counts)
  
  # Skip if insufficient groups
  if (length(group_counts) < 2) {
    warning("Insufficient groups in experiment ", exp)
    return(NULL)
  }
  
  # Step 1: Save the first column (HGNC symbols)
  hgnc_symbols <- residuals[[1]]  # or residuals_df$hgnc_symbol
  
  # Step 2: Subset columns that match row names of meta_exp
  matching_samples <- intersect(colnames(residuals), rownames(meta_exp))
  filtered_residuals <- residuals[, matching_samples, drop = FALSE]
  
  # Step 3: Combine HGNC symbols back with filtered residuals
  final_residuals <- cbind(hgnc_symbol = hgnc_symbols, filtered_residuals)
  
  # Define valid comparisons
  comparisons <- list(
    WT_vs_severe = c("WT", "severe_GoF"), # WT first
    WT_vs_mild = c("WT", "mild_GoF"), # WT first
    mild_vs_severe = c("mild_GoF", "severe_GoF") # mild first
  )
  
  # now we subset on comparisons
  for (comp_name in names(comparisons)) {
    groups <- comparisons[[comp_name]]
    
    meta_comp <- meta_exp %>%
      filter(group %in% groups) %>%
      mutate(group = factor(group, levels = groups))  # Relevel to 2 groups
    
    # quality check for appropriate group number
    if (n_distinct(meta_comp$group) < 2) {
      message(paste("Skipping", comp_name, "- insufficient groups"))
      next
    }
    
    # # Model fitting 
    # design <- model.matrix(~ group, data = meta_comp) ###### POTENTIAL PROBLEM AREA?
    # 
    # Proper design matrix
    design <- model.matrix(~ group, data = meta_comp) 
    cat("\nDesign matrix for", comp_name, ":\n")
    print(head(design))
    
    fit <- lmFit(final_residuals[, rownames(meta_comp)], design) %>%
      eBayes()
    print(design)
    
    # Get correct coefficient name
    coef_name <- colnames(design)[2]  # Will be "groupsevere_GoF" etc.
    cat("Testing coefficient:", coef_name, "\n")
    
    results <- tibble(
      gene_id = rownames(fit$coefficients),
      hgnc_symbol = hgnc_symbols,  # make sure this is in the same order!
      logFC = fit$coefficients[, coef_name],
      pvalue = fit$p.value[, coef_name],
      fdr = p.adjust(fit$p.value[, coef_name], method = "BH"),
      significant = fdr < fdr_threshold
    )
    
    results$experiment <- paste0("Exp", exp)
    results$comparison <- comp_name
    all_results[[paste0("Exp", exp, "_", comp_name)]] <- results
    
    # Store counts
    sig_count <- sum(results$significant, na.rm = TRUE)
    significant_counts[[paste0("Exp", exp, "_", comp_name)]] <- sig_count
    
    # Save results
    write.csv(
      results,
      file = paste0(results_dir, "DGE_Exp", exp, "_", comp_name, "_fdr_", fdr_threshold, ".csv"),
      row.names = FALSE)
    
    # Add metadata so we can identify which result this is from
    results$experiment <- exp
    results$comparison <- comp_name
    
    # Append to list
    all_results[[paste0("Exp", exp, "_", comp_name)]] <- results
    
    
    # Diagnostic MA plot ## chatgpt
    pdf(paste0(results_dir, "MA_", exp, "_", comp_name, ".pdf"))
    limma::plotMA(fit, main = paste(exp, comp_name))
    abline(h = 0, col = "red")
    dev.off()
    
  }
}

# Save significant counts ----
significant_counts_df <- enframe(significant_counts, name = "comparison", value = "n_sig") %>%
  separate(comparison, into = c("experiment", "contrast"), sep = "_", extra = "merge")


# SUMMARIZE RESULTS ----
# Convert significant counts to a tidy data frame
significant_counts_df <- significant_counts %>%
  enframe(name = "comparison", value = "significant_genes") %>%
  separate(comparison, into = c("experiment", "contrast"), sep = "_", extra = "merge") %>%
  pivot_wider(names_from = contrast, values_from = significant_genes)

print(significant_counts_df)

# Combine all results
final_results <- bind_rows(all_results)

# Save results
write.csv(final_results, paste0(results_dir, 
                                "DGE_results_all_comparisons", 
                                "_fdr_", fdr_threshold, 
                                ".csv"), row.names = FALSE)
########

# Reshape to long format
significant_counts_df <- pivot_longer(significant_counts_df, 
                                      cols = -experiment, 
                                      names_to = "contrast", 
                                      values_to = "number_significant")

# Optional: set factor levels for prettier facet ordering
significant_counts_df$contrast <- factor(significant_counts_df$contrast, levels = c("WT_vs_severe", "WT_vs_mild", "mild_vs_severe"))

# Ensure numeric y-axis
significant_counts_df$number_significant <- as.numeric(significant_counts_df$number_significant)

pdf("number_of_DEgenes_by_experiment.pdf", width = 5, height = 3.5)

significant_counts_df$number_significant <- unlist(significant_counts_df$number_significant)

ggplot(significant_counts_df, aes(x = experiment, y = number_significant, fill = experiment)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_text(aes(label = number_significant), 
            vjust = 0, size = 3, position = position_dodge(0.9)) + 
  facet_wrap(~ contrast, nrow = 1) + 
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_blank()
  ) + 
  ylab("Number of DE Genes") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Corrected to continuous scale
  labs(fill = "Experiment")

dev.off()

de_gene_table <- results %>%
  filter(significant) %>%
  select(hgnc_symbol, logFC, pvalue, fdr)

write.csv(de_gene_table, "DE_genes_with_values.csv", row.names = FALSE)

# Filter for only significant DE genes
de_genes_all <- results %>%
  filter(significant) %>%
  select(hgnc_symbol, logFC, fdr, experiment, comparison)

# Optional: clean up labels for plot
de_genes_all$comparison <- factor(de_genes_all$comparison,
                                  levels = c("WT_vs_severe", "WT_vs_mild", "mild_vs_severe"),
                                  labels = c("WT vs Severe", "WT vs Mild", "Mild vs Severe"))

library(dplyr)
library(ggplot2)

# Ensure experiment is a factor
final_results$experiment <- as.factor(final_results$experiment)

# Count significant DE genes per experiment and comparison
sig_counts <- final_results %>%
  filter(significant == TRUE) %>%
  group_by(experiment, comparison) 

# Save results
write.csv(sig_counts, paste0(results_dir, 
                             "DEgenes_by_experiment_comparisons", 
                             "_fdr_", fdr_threshold, 
                             ".csv"), row.names = FALSE)

library(ggplot2)
library(dplyr)

# Combine experiment and comparison into a single label for plotting
sig_counts$label <- paste0("Exp ", sig_counts$experiment, ": ", sig_counts$comparison)

# Ensure label is a factor for consistent facet order
sig_counts$label <- factor(sig_counts$label)

library(ggplot2)
library(dplyr)

# Create faceting label
sig_counts <- sig_counts %>%
  mutate(
    experiment = as.factor(experiment),
    comparison = as.factor(comparison),
    label = paste0("Exp ", experiment, ": ", comparison),
    significant = as.logical(significant)
  )

# Plot
p <- ggplot(sig_counts, aes(x = logFC, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  facet_wrap(~label, scales = "free") +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey70")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Volcano Plot by Experiment and Comparison",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significant"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
pdf("test_volcano.pdf", width = 10, height = 7)
print(p)
dev.off()