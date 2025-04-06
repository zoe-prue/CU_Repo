## linux-64
salloc
module load conda

## create and activate environment
conda create -n mapping_env ## only once initially

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

## check versions
fastqc -v
trim_galore --version
kallisto version
multiqc --version

## just run fastqc on raw fastqs
## need to call fastqc
## haley provides script edited to be comptaible with my data
