cd /groups/db3700_gp/zp2280/raw_data/RNAseq_Exp21-25-26/30-1116451817/00_fastq

## say i didn't know experimental design, i can count the number of outputs by # of lines
ls | grep fastq.gz | wc -l


## one-linder script = running a for loop in a single line
## this one linder will generate and submit a script for each fastq
## this one is just for the first fastq
## submit into terminal
for i in *.fastq.gz;do echo -e '#!/bin/sh\nsource activate base;conda activate mapping_env; fastqc' $i>$i.sh;sbatch --time=1:0:0 --mem=4G --nodes=1 $i.sh;done
## time -> much longer than job would take, but not too long 
	# bc it would be pushed back due to lack of resources immediately available
## creates a shell script named i.sh
## one liner is semicolons

## check that job is running (make sure ST is R)
squeue -u zp2280

## less to open file
	# see the i.sh script could be helpful
## outputs for fastqc are .hmtl, .zip, i.sh, and slurm.ou
## instead of looking at every single output by hand, make it a wildcard
## looking for "Analysis complete" 
## last ten lines, whats the word count (wc) and line number
tail *out | grep 'Analysis complete' | wc -l # shows last ten lines (54)

## can cancel a job with 
scancel

## making a folder for each type of output
mkdir scripts
mkdir files_for_multiqc
mkdir slurm_outs

## moving outputs
mv *.sh scripts/
mv *fastqc.zip files_for_multiqc
mv *fastqc.html files_for_multiqc
mv *.out slurm_outs

## go into the files_for_multiqc
	# aggrgate the fastqc outouts into multiqc
	multiqc .
	mv multiqc_report.html multiqc_pretrimming.html

## go to Globus, transfer this file to my computer and open (it's a webpage)
	# how many reads per sample do i have
	# do you have round the same number of reads per sample
	# Q30 quality of bases
	# look up phred scores = read quality metric

## do not run, but know the concepts
#!/bin/sh # run before every shell script
# no \n but that just means new line
source activate base
conda activate mapping_env
fastqc xx.fasq.gz # fastq file name individually

## this is to read easily; not a one liner bc no semicolon
for i in Ep26-JAK1-E139K-1_R1_001.fastq.gz
	do echo -e '#!/bin/sh\nsource activate base;conda activate mapping_env; fastqc' $i>$i.sh
	sbatch --time=1:0:0 --mem=4G --nodes=1 $i.sh
done



