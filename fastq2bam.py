#! /bin/env python
# npagane | 190618 | risca lab | execute fastq2bam pipeline remotely

import os
import sys
import time

start = time.time()

# parameters needed to insert as argruments to execute
inputs = [
    # how many cores to run (integer, i.e. 14)
    'core',
    # directory with fastq files (string, i.e. path/to/fastq)
    'fastqDir'
]

# optional tags
optTags = {
    # genome reference for alignment (string, i.e. path/to/genome.fa)
    '--genomeRef': "/rugpfs/fs0/risc_lab/scratch/nvelez/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome",
    # blacklist for filtering (string, i.e. path/to/blacklist)
    '--blacklist': "/rugpfs/fs0/risc_lab/scratch/nvelez/blacklists/ATAC_blacklist.bed",
    # the set number of the experiment (int, i.e. 001)
    '--set': '001', 
    # the map quality threshold for alignment (int, i.e. 30)
    '--mapq': '30', 
    # whether there are fastq files for the index reads or not
    '--index': 'True',
    # any snakemake flags for compilation (string, i.e. unlock)
    '--snakemake': '', 
   
}

if __name__ == '__main__':
    # make sure user has specified number of cores and working directory
    if len(sys.argv)-1 < 2:
       print('need at least 2 arguments: how many cores to run pipeline and where the fastq files are (i.e. working directory)')
       sys.exit(0)
    # gather number of cores and working directory
    cores = sys.argv[1]
    wd = sys.argv[2]
    # look for optional tags
    addedTags = ''
    for tag in optTags.keys():
        tempTag = ' '.join(sys.argv).split(tag)
        if len(tempTag) == 1:
            addedTags += " '" + tag[2:] + '="' + optTags[tag] + '"' + "'"
        else:
            tempFlag = tempTag[1].split(' ')[1]
            if tag == "--snakemake":
                addedTags += " --" + tempFlag
            else:
                addedTags += " '" + tag[2:] + '="' + tempFlag + '"' + "'"
    # set working directory to where the Snakefile is
    os.chdir('/rugpfs/fs0/risc_lab/store/npagane/fastq2bam')
    # execute snakemake in actual working directory and with proper configs
    os.system('snakemake --rerun-incomplete --cores ' + cores + ' --directory ' + wd + ' --config' + addedTags)
    stop = time.time()
    print('ran took ' + str(1.0*(stop - start)/(60*60)) + ' hours')
