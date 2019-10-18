#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --array=1-1000


SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p sample.txt)
echo $SAMPLE

# run code
python ~/exe/fastq2bam.py 2 . /rugpfs/fs0/risc_lab/store/risc_data/raw-data/RICC-testing/AndrewScortea/MiSeq_07122019 $SAMPLE

python 

