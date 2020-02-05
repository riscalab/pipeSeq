#!/bin/bash
# npagane | risca lab | oct 2019 | ATACseq pipeline wrapper to execute

##################
# SET PARAMETERS #
##################

# set the default arguments for optional parameters
genomeMap="hg38"
mapq=30
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
g) genomeMap=${OPTARG};; # genome reference for alignment (OPTIONAL, string, OPTS: 'hg38', 'hg19', 'mm9', 'mm10')
b) blacklist=${OPTARG};; # blacklist for filtering (OPTIONAL, string, defaults to genomeMap blacklist OR overwrite with path/to/blacklist)
m) mapq=${OPTARG};; # the map quality threshold for alignment (OPTIONAL, int, i.e. 30)
p) snakemake=${OPTARG};; # any snakemake flags for compilation (OPTIONAL, string, i.e. "--snakemake unlock")
esac
done

# check for required arguments
if [ -z "$cwd" ] || [ -z "$fastqDir" ] || [ -z "$sampleText" ]
then
    echo "must specificy the desired working directory (c), fastq directory (f), and text file with sample names (s)"
    exit
fi

# check genomeMap to align to correct genome with either assumed or specified blacklist
if [ "$genomeMap" != "hg38" ]
then 
    if [ "$genomeMap"== "hg19" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome"
        TSS="/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/hg19/GENCODE/gencode.v19.tss.bed"
        chromSize="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/chrom.sizes"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/ATAC_blacklist.bed"
        fi
    elif  [ "$genomeMap"== "mm9" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/genome/Sequence/Bowtie2Index/genome"
        TSS="/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/mm9/GENCODE/mouse.gencode.m4.tss.bed"
        chromSize="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/genome/chrom.sizes"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/blacklist/ATAC_blacklist.bed"
        fi
    elif [ "$genomeMap"== "mm10" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/Bowtie2Index/genome"
        #TSS=""
        #chromSize=""
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/blacklist/ATAC_blacklist.bed"
        fi
    fi
else
    genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome"
    TSS="/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/hg38/GENCODE/gencode.v30.basic.annotation.tss.bed"
    chromSize="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes"
    if [ -z "$blacklist" ]
    then
        blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/ATAC_blacklist.bed"
    fi
fi
echo "alinging to $genomeRef with blacklist $blacklist with chromosome sizes $chromSize and TSS locations $TSS"

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
