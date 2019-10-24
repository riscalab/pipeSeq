#!/bin/bash
# npagane | risca lab | oct 2019 | wrapper to unlock all the relevant snakemake directories

# run this script just like you would run the ATACseq.sh exectubale
# i.e. bash unlock.sh -c . -f /path/to/fastqs -s samples.txt

##################
# SET PARAMETERS #
##################

# set the default arguments for optional parameters
genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome"
blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/ATAC_blacklist.bed"
TSS="/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/hg38/GENCODE/gencode.v30.basic.annotation.tss.bed"
mapq=30
chromSize="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes"
snakemake="--unlock"
# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/ATACseq"
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/ATACseq"

# parse the arguments
while getopts c:f:s:g:b:t:m:z:p: option
do
case "${option}"
in
c) cwd=${OPTARG};; # working directory for analysis (REQUIRED, i.e. path/to/workingDirectory)
f) fastqDir=${OPTARG};; # directory with fastq files (REQUIRED, string, i.e. path/to/fastq)
s) sampleText=${OPTARG};; # sample names text file (REQUIRED, string, i.e. path/to/samples.txt)
g) genomeRef=${OPTARG};; # genome reference for alignment (OPTIONAL, string, i.e. path/to/genome.fa)
b) blacklist=${OPTARG};; # blacklist for filtering (OPTIONAL, string, i.e. path/to/blacklist)
t) TSS=${OPTARG};; # TSS for that genome (OPTIONAL, string, i.e. path/to/TSS) 
m) mapq=${OPTARG};; # the map quality threshold for alignment (OPTIONAL, int, i.e. 30)
z) chromSize=${OPTARG};;  # sizes of chromosomes (string, i.e. path/to/genome/file)
p) snakemake=${OPTARG};; # any snakemake flags for compilation (OPTIONAL, string, i.e. "--snakemake unlock")
esac
done

# check for required arguments
if [ -z "$cwd" ] || [ -z "$fastqDir" ] || [ -z "$sampleText" ]
then
    echo "must specificy the desired working directory (c), fastq directory (f), and text file with sample names (s)"
    exit
fi

source activate ATACseq

# unlock directories
for i in `cat $sampleText`
do
    snakemake --snakefile $exeDir/Snakefile --rerun-incomplete --cores 1 $snakemake --config "fastqDir='$fastqDir'" "genomeRef='$genomeRef'" "blacklist='$blacklist'" "TSS='$TSS'" "mapq='$mapq'" "chromSize='$chromSize'" "sample='$i'"
done
