#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4

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
singleend=$7
exeDir=$8
snakemake=$9

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
snakemake --snakefile $exeDir/Snakefile --nolock --rerun-incomplete --cores 4 $snakemake --config "pipe='$pipe'" "fastqDir='$fastqDir'" "genomeRef='$genomeRef'" "blacklist='$blacklist'" "mapq='$mapq'" "sample='$SAMPLE'" "singleend='$singleend'"
