### BEFORE: putting star into mapping_env
## linux-64
salloc
module load conda

## next time, dont "install" mapping_env again, just do this
source activate base ## always need to run to activate
conda activate mapping_env ## always need to run to activate

# conda install -c bioconda (in order)
conda install -c bioconda -n mapping_env trim-galore
conda install -c bioconda -n mapping_env fastqc
conda install -c bioconda -n mapping_env multiqc
conda install -c bioconda -n mapping_env kallisto
conda install -c bioconda -n mapping_env STAR
conda list -n mapping_env



#!/bin/bash
## to generate scripts for star
cd /groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/trimmed_fastqs
OUTS_DIR=/groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/aligned_fastqs

for FASTQ in *_R1_001_val_1.fq.gz
do
	SAMPLE=${FASTQ%%_*}

	echo -e "#!/bin/bash\nmodule load conda\nsource activate base\nconda activate mapping_env\nSTAR --runThreadN 8 --genomeDir /groups/db3700_gp/SHARED/STAR_indexes/genomeDir/ --readFilesIn "${SAMPLE}"_R1_001_val_1.fq.gz "${SAMPLE}"_R2_001_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix "${OUTS_DIR}"/star_PE_trimgalore_PE_"${SAMPLE}" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts" > ${SAMPLE}.star_PE_trimgalore_PE.sh

	sbatch --time=01:00:00 --nodes=1 --cpus-per-task=8 --mem-per-cpu=32G ${SAMPLE}.star_PE_trimgalore_PE.sh
done


### index dirctory and file: /groups/db3700_gp/SHARED/kallisto_indexes/Homo_sapiens.GRCh38.100.v13.idx

##/groups/db3700_gp/SHARED/STAR_indexes/genomeDir/ --readFilesIn /groups/db3700_gp/er3277/HD1/FASTQ/TRIMMED/61-1_R1_001_val_1.fq.gz groups/db3700_gp/er3277/HD1/FASTQ/TRIMMED/61-1_R2_001_val_2.fq.gz --readFilesCommand zcat --outFileNamePrefix /groups/db3700_gp/er3277/HD1/STARout/1/ --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

#/groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/aligned_fastqs/ star_PE_trimgalore_PE_Exp21-JAK1-S703I-2_R1_001_val_1.fq.gz
#/star_PE_trimgalore_PE_"${FASTQ}