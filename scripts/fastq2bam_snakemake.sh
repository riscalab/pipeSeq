#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3

######################
# READ IN PARAMETERS #
######################

pipe="fastq2bam"
cwd=$1
fastqDir=$2
sampleText=$3
genomeRef=$4
blacklist=$5
mapq=$6
snakemake=$7
# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/pipeSeq"
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/pipeSeq"

###########
# EXECUTE #
###########

echo This is job $SLURM_JOB_ID

# conda env
echo `which python`

# determine sample
SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleText)
echo $SAMPLE

# run code
snakemake --snakefile $exeDir/Snakefile --nolock --rerun-incomplete --cores 3 $snakemake --config "pipe='$pipe'" "fastqDir='$fastqDir'" "genomeRef='$genomeRef'" "blacklist='$blacklist'" "mapq='$mapq'" "sample='$SAMPLE'" 

