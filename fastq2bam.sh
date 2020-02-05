#!/bin/bash
# npagane | risca lab | oct 2019 | fastq2bam pipeline wrapper to execute

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
while getopts c:f:s:g:b:m:p: option
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
    if [ "$genomeMap" == "hg19" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed"
        fi
    elif  [ "$genomeMap" == "mm9" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/blacklist/mm9-blacklist.bed"
        fi
    elif [ "$genomeMap" == "mm10" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/blacklist/mm10-blacklist.v2.bed"
        fi
    else
        echo "unrecognized genome.\navailable genomes: hg38, hg19, mm10, mm9.\ntalk to nicole to get your genome on the cluster if not there.\n"
        exit
    fi
else
    genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome"
    if [ -z "$blacklist" ]
    then
        blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/hg38-blacklist.v2.bed"
    fi
fi
echo "alinging to $genomeRef with blacklist $blacklist"

# change working directory
cd $cwd

# set conda env
source activate fastq2bam

# run fastq2bam
numSamples=`wc -l $sampleText | awk '{print $1}' `
fastq2bam=$(sbatch -p risc,hpc --array=1-$numSamples $exeDir/scripts/fastq2bam_snakemake.sh $cwd $fastqDir $sampleText $genomeRef $blacklist $mapq $snakemake)

# get job id
if ! echo ${fastq2bam} | grep -q "[1-9][0-9]*$"; then
   echo "Job(s) submission failed."
   echo ${fastq2bam}
   exit 1
else
   fastq2bamID=$(echo ${fastq2bam} | grep -oh "[1-9][0-9]*$")
fi

# summary stats for fastq2bam after successful completion
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
sbatch -p risc,hpc --depend=afterok:$fastq2bamID --wrap="python $exeDir/rules/helper.py 0 $fastq2bamID $sampleText '$myInvocation' $fastqDir $genomeRef $blacklist $mapq"
