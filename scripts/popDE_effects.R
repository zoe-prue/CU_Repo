library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggfortify)
ul = function(x){x[1:6, 1:6]}

setwd("/Users/Haley/Desktop/lab/code/MAZ_FLU/")
working_dir <- "effect_of_KD/"
model_folder <- "popDE/"
dir.create(paste0(working_dir,model_folder))
name_new_folder <- "results/"
dir.create(paste0(working_dir,model_folder,name_new_folder))

## outputs
results_dir <- paste0(working_dir,model_folder,name_new_folder)

## read in meta data
meta_data <- read.csv("meta_data.csv", header = TRUE)
meta_data$infection_status = factor(meta_data$infection_status, levels = c("NI","flu"))
meta_data$infection_mouse <- paste0(meta_data$mouse_type,"_",meta_data$infection_status)
# factor WT, then mild, then severe
meta_data$infection_mouse = factor(meta_data$infection_mouse, levels = c("WT", "WT_flu", "CypDko_NI", "CypDko_flu"))
meta_data$mouse_type = factor(meta_data$mouse_type, levels = c("WT","CypDko"))
rownames(meta_data) <- meta_data$sample_ID

## read in expression data, not batch corrected with without JAK1-WT_exp26_3
# log2cpm with batch effect included
file <- "GE_matrix_corrected_MalignedScale_logCPM_voomed"
weights_file <- "GE_matrix_corrected_MalignedScale_weights"

residuals <- read.table(paste0("log2cpm_expression/",file,".txt"), header = TRUE, sep = ",")
weights <- read.table(paste0("log2cpm_expression/",weights_file,".txt"), header = TRUE, sep = ",")

length(which(colnames(residuals)!=rownames(meta_data)))
length(which(colnames(weights)!=rownames(meta_data)))
length(which(colnames(residuals)!=colnames(weights)))

## subset the meta data on only the WT and severe samples
## refactor the Group to only those 2 levels
## subset the expression data on those same samples
## make sure the rownames of the meta data match the colnames of the expression data 

## fit the model
design = model.matrix(~ batch + infection_status, data = meta_data)
vfit <- lmFit(residuals, weights = weights, design)
vfit <- eBayes(vfit)

head(vfit$coefficients)

## get results
betas = as.data.frame(vfit$coefficients[, which(colnames(vfit$coefficients) %in% c("GroupKO"))]); colnames(betas)[1] <- "betas"
p_values = as.data.frame(vfit$p.value[, which(colnames(vfit$coefficients) %in% c("GroupKO"))]); colnames(p_values)[1] <- "pvalues"
fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH")); colnames(fdrs)[1] <- "fdrs"

results <- cbind(betas, p_values, fdrs)
write.table(results, paste0(results_dir,"resultsALL.txt"), quote = FALSE, sep = ",")

## get significant gene numbers
numDEgenes_20 <- length(which(fdrs[,1] < 0.20))
numDEgenes_15 <- length(which(fdrs[,1] < 0.15))
numDEgenes_10 <- length(which(fdrs[,1] < 0.10))
numDEgenes_05 <- length(which(fdrs[,1] < 0.05))

numDEgenes <- as.data.frame(t(as.data.frame(cbind(numDEgenes_20, numDEgenes_15, numDEgenes_10, numDEgenes_05, numDEgenes_20_flu, numDEgenes_15_flu, numDEgenes_10_flu, numDEgenes_05_flu))))
colnames(numDEgenes)[1] <- "number_significant"
write.table(numDEgenes, paste0(results_dir,"num_popDEgenes.txt"), quote = FALSE, sep = ",")

## DE genes
## to plot
numDEgenes$variable <- rownames(numDEgenes)
numDEgenes_p <- numDEgenes

## parse infection status
numDEgenes_p <- separate(numDEgenes_p, variable, sep = "_", remove = FALSE, into = c("x","FDR_cutoff","infection"))
numDEgenes_p$number_significant <- as.numeric(numDEgenes_p$number_significant)
numDEgenes_p$FDR_cutoff <- as.numeric(numDEgenes_p$FDR_cutoff)
numDEgenes_p$FDR_cutoff <- (numDEgenes_p$FDR_cutoff)/100
numDEgenes_p$FDR_cutoff <- as.factor(numDEgenes_p$FDR_cutoff)
numDEgenes_p$infection = factor(numDEgenes_p$infection, levels=c("NI","flu"))

## plots
pdf(paste0(results_dir,"number_of_popDEgenes.pdf"), width = 5, height = 3.5)
ggplot(numDEgenes_p, aes(fill = FDR_cutoff, y = number_significant, x = x)) +
	geom_bar(position="dodge", stat="identity") + 
	geom_text(aes(label = number_significant, color = FDR_cutoff), position = position_dodge(1), vjust = 0, size = 2.5) +
	facet_wrap(vars(infection), nrow = 1) +
	#theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
		  panel.grid.major = element_blank(),
		  panel.background = element_rect(fill = "white"),
		  strip.background = element_rect(fill = NA)) +
	scale_fill_manual(values = c("#3D3242", "#7A3E48", "#B95835", "#E18942"), labels = c("0.05","0.10","0.15","0.20")) +
	scale_color_manual(values = c("#3D3242", "#7A3E48", "#B95835", "#E18942")) + 
	labs(fill = "FDR") +
	guides(color = FALSE) +
	xlab("") +
	scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
dev.off()


