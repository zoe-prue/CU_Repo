cd /groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/trimmed_fastqs

for i in *.fq.gz;do echo -e '#!/bin/sh\nmodule load conda\nsource activate base;conda activate mapping_env;fastqc' $i>$i.sh;sbatch --time=1:00:0 --mem=4G --nodes=1 $i.sh;done

## making a folder for each type of output
mkdir scripts
mkdir files_for_multiqc_post_trims
mkdir slurm_outs_post_trim

## moving outputs
mv *.sh scripts/
mv *fastqc.zip files_for_multiqc_post_trim
mv *fastqc.html files_for_multiqc_post_trim
mv *.out slurm_outs_post_trim
mv *.fastq.gz_trimming_report.txt files_for_multiqc_post_trim

## go into the files_for_multiqc
	# aggrgate the fastqc outouts into multiqc
	multiqc .
	mv multiqc_report_1.html multiqc_posttrimming.html

# 15:45:24 zp2280@login1:files_for_multiqc_post_trim# module load conda
# 15:45:49 zp2280@login1:files_for_multiqc_post_trim# source activate base
# (base) 15:45:59 zp2280@login1:files_for_multiqc_post_trim# conda activate mapping_env
# (mapping_env) 15:46:13 zp2280@login1:files_for_multiqc_post_trim# multiqc .

