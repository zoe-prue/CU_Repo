#!/bin/bash
## to generate scripts for kallisto
cd /groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/trimmed_fastqs
OUTS_DIR=/groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/aligned_fastqs

for FASTQ in Ep26-JAK1-E139K-1_R1_001_val_1.fq.gz
do
	FASTQ_NAME=${FASTQ/\_R*/}

	echo -n "#!/bin/bash
	
	module load conda
	source activate base
	conda activate mapping_env

	kallisto quant -i /groups/db3700_gp/SHARED/kallisto_indexes/Homo_sapiens.GRCh38.100.v13.idx -o "${OUTS_DIR}"/kallisto_PE_trimgalore_PE_"${FASTQ}" "${FASTQ}"_L001_R1_001_val_1.fq.gz "${FASTQ}"_L001_R2_001_val_2.fq.gz" > ${FASTQ}.kallisto_PE_trimgalore_PE.sh

	sbatch --time=5:00:0 --mem=20G --nodes=1 ${FASTQ}.kallisto_PE_trimgalore_PE.sh
done


### index dirctory and file: /groups/db3700_gp/SHARED/kallisto_indexes/Homo_sapiens.GRCh38.100.v13.idx