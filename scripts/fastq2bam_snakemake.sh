#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3

######################
# READ IN PARAMETERS #
######################

cwd=$1
fastqDir=$2
sampleText=$3
genomeRef=$4
blacklist=$5
TSS=$6
mapq=$7
snakemake=$8
# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/fastq2bam"
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/fastq2bam"

###########
# EXECUTE #
###########

echo This is job $SLURM_JOB_ID
echo Test permission

# conda env
source activate fastq2bam 
echo `which python`

# determine sample
SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleText)
echo $SAMPLE

# run code
snakemake --snakefile $exeDir/Snakefile --rerun-incomplete --cores 3 $snakemake --config "fastqDir='$fastqDir'" "genomeRef='$genomeRef'" "blacklist='$blacklist'" "TSS='$TSS'" "mapq='$mapq'" "sample='$SAMPLE'" 

