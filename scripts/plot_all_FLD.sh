#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p risc,hpc

# get in the correct conda env
source activate fastq2bam

# read in libs
libs=$1

# make temp FLD to read in 
for i in `cat $libs`
do
   samtools view $i/*rmdup.bam | awk '{print $9}' > $i/"tempFLD.log"
done

# run FLD script
Rscript /rugpfs/fs0/risc_lab/store/npagane/pipeSeq/scripts/plotisds_v3.R $libs tempFLD

# remove temp FLD 
for i in `cat $libs`
do
   rm $i/"tempFLD.log"
done
