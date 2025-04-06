library(biomaRt)
library(plyr)
library(dplyr)
ul <- function(df){df[1:6,1:6]}

setwd("~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/counts_kallisto")

## in excel, split into 3 experiments before reading in 3 .txt?

## read in counts data
GE <- read.table(paste0("Homo_sapiens.GRCh38.87.kallisto.gene_level.lengthScaledTPM_counts.txt"), header = TRUE)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
# Alternative mirrors:
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# ?useEnsembl()

colnames(GE)  # To see a sample of the row names

genes_dups_ensembl <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"),    
                      filters = "ensembl_gene_id",
                      values = rownames(GE),
                      mart = ensembl)

genes_dups_ensembl <- unique(genes_dups_ensembl)

# table(genes_dups_ensembl$ensembl_gene_id)

## only keep unique mgi ids
genes_ensembl <- genes_dups_ensembl[!duplicated(genes_dups_ensembl$hgnc_symbol), ]

## add as col
GE$ensembl_gene_id <- rownames(GE)

## merge unique genes and dataset
GE <- join(genes_ensembl, GE, by = "ensembl_gene_id", type = "left")
## 44787 genes before low expression filtering

## Number of genes that are missing an external_gene_name
which(is.na(GE$mgi_symbol)) %>% length

rownames(GE) <- GE$mgi_symbol
GE$ensembl_gene_id <- NULL
GE$mgi_symbol <- NULL

write.table(GE, "Homo_sapiens.GRCh38.87.kallisto.gene_level.lengthScaledTPM_counts_HGNCsymbols.txt", quote = FALSE, sep = ",")


