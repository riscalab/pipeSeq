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
chromSize=$8
snakemake=$9
# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/ATACseq" #"/rugpfs/fs0/risc_lab/store/npagane/ATACseq" 

###########
# EXECUTE #
###########

# conda env
source activate ATACseq
echo `which python`

# determine sample
SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sampleText)
echo $SAMPLE

# run code
cd $cwd
snakemake --snakefile $exeDir/Snakefile --rerun-incomplete --cores 3 $snakemake --config "fastqDir='$fastqDir'" "genomeRef='$genomeRef'" "blacklist='$blacklist'" "TSS='$TSS'" "mapq='$mapq'" "chromSize='$chromSize'" "sample='$SAMPLE'" 

