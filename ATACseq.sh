#!/bin/bash
# npagane | risca lab | oct 2019 | ATACseq pipeline wrapper to execute

##################
# SET PARAMETERS #
##################

# set the default arguments for optional parameters
genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome"
blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/ATAC_blacklist.bed"
TSS="/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/hg38/GENCODE/gencode.v30.basic.annotation.tss.bed"
mapq=30
chromSize="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes"
snakemake=""
# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/pipeSeq" 
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/pipeSeq"

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

# change working directory
cd $cwd

# set conda env
source activate ATACseq

# run ATACseq
numSamples=`wc -l $sampleText | awk '{print $1}' `
ATACseq=$(sbatch -p risc,hpc --array=1-$numSamples $exeDir/scripts/ATACseq_snakemake.sh $cwd $fastqDir $sampleText $genomeRef $blacklist $TSS $mapq $chromSize $snakemake)

# get job id
if ! echo ${ATACseq} | grep -q "[1-9][0-9]*$"; then
   echo "Job(s) submission failed."
   echo ${ATACseq}
   exit 1
else
   ATACseqID=$(echo ${ATACseq} | grep -oh "[1-9][0-9]*$")
fi

# summary stats for fastq2bam after successful completion
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
sbatch -p risc,hpc --depend=afterok:$ATACseqID --wrap="python $exeDir/rules/helper.py 1 $ATACseqID $sampleText '$myInvocation' $fastqDir $genomeRef $blacklist $mapq $TSS $chromSize"
