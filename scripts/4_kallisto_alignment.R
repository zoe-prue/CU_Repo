# count matrices -> convert ENS's to hgncs -> removal of lowly expressed genes using limma voom -> PCAs -> DGE

# module load R
# R

library(GenomicFeatures)
library(tximport)
library(rhdf5)

# Informations for output file names 
organism='Homo_sapiens'
version='GRCh38.87'
aligner='kallisto'
OUT_DIR='~/Desktop/CU_coding/RNA-seq/DGE_RNAseq/counts_kallisto'

# Read sample informations/quant files location 
gtf_filename="~/Desktop/CU_coding/RNA-seq/meta_data/Homo_sapiens.GRCh38.87.gtf"
#sample_filename=""

samples=list.files("~/Desktop/CU_coding/RNA-seq/raw_kallisto_alignment")
samples <- gsub("kallisto_PE_trimgalore_PE_","",samples)
files=samples
files=paste0('~/Desktop/CU_coding/RNA-seq/raw_kallisto_alignment/kallisto_PE_trimgalore_PE_',files,'/abundance.h5')
names(files)=samples

# Build GeneDB
# could give me issues due to version controls
txdb=makeTxDbFromGFF(gtf_filename,format='gtf')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Transcript-level raw counts
txi.kallisto.transcript_level.counts <- tximport(files, type = "kallisto", txOut = TRUE)

# Transcript-level scaledTPM
txi.kallisto.transcript_level.scaledTPM <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance= 'scaledTPM')

# Transcript-level lengthScaledTPM
txi.kallisto.transcript_level.lengthScaledTPM <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance= 'lengthScaledTPM') # scaled by transcript length; what haley uses

# Gene-level raw counts
txi.kallisto.gene_level.counts <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene, ignoreTxVersion = TRUE)

# Gene-level scaledTPM
txi.kallisto.gene_level.scaledTPM <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene, countsFromAbundance= 'scaledTPM', ignoreTxVersion = TRUE)

# Gene-level lengthScaledTPM
txi.kallisto.gene_level.lengthScaledTPM <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene, countsFromAbundance= 'lengthScaledTPM', ignoreTxVersion = TRUE)


# Write tables 
write.table(txi.kallisto.transcript_level.counts, paste0(OUT_DIR,"/",paste(organism,version,aligner,'transcript_level','counts','txt',sep='.')),quote=F)
write.table(txi.kallisto.transcript_level.scaledTPM, paste0(OUT_DIR,"/",paste(organism,version,aligner,'transcript_level','scaledTPM','txt',sep='.')),quote=F)
write.table(txi.kallisto.transcript_level.lengthScaledTPM, paste0(OUT_DIR,"/",paste(organism,version,aligner,'transcript_level','lengthScaledTPM','txt',sep='.')),quote=F)

write.table(txi.kallisto.gene_level.counts, paste0(OUT_DIR,"/",paste(organism,version,aligner,'gene_level','counts','txt',sep='.')),quote=F)
write.table(txi.kallisto.gene_level.scaledTPM, paste0(OUT_DIR,"/",paste(organism,version,aligner,'gene_level','scaledTPM','txt',sep='.')),quote=F)
write.table(txi.kallisto.gene_level.lengthScaledTPM, paste0(OUT_DIR,"/",paste(organism,version,aligner,'gene_level','lengthScaledTPM','txt',sep='.')),quote=F)

# the ones i want with counts
write.table(txi.kallisto.transcript_level.counts$counts, paste0(OUT_DIR,"/",paste(organism,version,aligner,'gene_level','raw_counts','txt',sep='.')),quote=F)

# generally used for PCA
write.table(txi.kallisto.gene_level.lengthScaledTPM$counts, paste0(OUT_DIR,"/",paste(organism,version,aligner,'gene_level','lengthScaledTPM_counts','txt',sep='.')),quote=F)




