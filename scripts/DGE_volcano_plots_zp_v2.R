library(ggplot2)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(data.table)
library(ggrepel)

setwd("/Users/Haley/Desktop/lab/code/COVID_19_MONTREAL_COMBINED")

results_dir <- paste0("results/231016_modeling/infection_effects/")
celltypes <- c("B","CD4_T","CD8_T","CD14_monocytes","CD16_monocytes","NK")

fdr_thresh = 0.05
beta_thresh = 0.5
plots_list <- list()

label <- function(values, beta_down, beta_up) {
  return(ifelse(values <= beta_down, 'downregulated', ifelse(values >= beta_up, 'upregulated', 'no change')))
}

for(i in (1:length(celltypes))){
  
  celltype_i <- celltypes[i]
  
  DE_results <- read.table(paste0(results_dir,celltype_i,"_results_with_qvalues.txt"), header = TRUE, sep = ",")
  DE_results$label <- label(DE_results[,"betas"], -0.5, 0.5)
  DE_results$label <- as.factor(DE_results$label)
  DE_results$genes <- rownames(DE_results)
  
  plots_list[[i]] <- ggplot(DE_results, aes(betas, -log10(pvalues))) +
    geom_point(color = "grey50", size = 2, alpha = 0.2) +
    geom_point(data = subset(DE_results, Fdr_SAB < fdr_thresh), aes(color = factor(label, levels = c("upregulated","no change","downregulated"))), alpha = 0.3) +
    geom_hline(yintercept= -log10(0.05), linetype="dashed", color="grey") +
    geom_vline(xintercept= c(-1 * beta_thresh, beta_thresh), linetype = "dashed", color = "grey60") +
    scale_colour_manual(values=c("#cd7058", "grey50", "#599ad3")) +
    xlab("log2FC") +
    ylab("-log10(p-value)") +
    ggtitle(paste0(celltype_i)) +
    theme_bw() +
    theme(plot.title = element_text(size = 15), legend.title=element_blank(),
          legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    geom_text_repel(data = subset(DE_results, 
                                  if(celltype_i == "CD14_monocytes" || celltype_i == "NK" || celltype_i == "CD16_monocytes"){
                                    betas > 2
                                  }else if(celltype_i == "B" || celltype_i == "CD8_T" || celltype_i == "CD4_T"){
                                    betas > 1
                                  }else{
                                    betas > 1
                                  }), 
                    aes(betas, -log10(pvalues), label = genes), color = "black", size = 3, hjust = 0, segment.size = 0.2, segment.color = "black", segment.alpha = 0.3, fontface = "italic") +
    geom_text_repel(data = subset(DE_results, 
                                  if(celltype_i == "CD14_monocytes" || celltype_i == "monocytes_combined" || celltype_i == "CD16_monocytes" || celltype_i == "NK"){
                                    betas < -1.75
                                  }else{
                                    betas < -1
                                  }), 
                    aes(betas, -log10(pvalues), label = genes), color = "black", size = 3, hjust = 0, segment.size = 0.2, segment.color = "black", segment.alpha = 0.3, fontface = "italic")
  
  print(celltype_i)
}

pdf(paste0(results_dir,"volcano_COVID_effect_FDR0.05_logFC0.5.pdf"), width = 14, height = 9)
grid.arrange(grobs = plots_list, ncol = 3)
dev.off()