#! /bin/env python
# npagane | 190618 | risca lab | execute ATACseq pipeline remotely

import os
import sys
import time

start = time.time()

# parameters needed to insert as argruments to execute
inputs = [
    # how many cores to run (integer, i.e. 14)
    'core',
    # directory with fastq files (string, i.e. path/to/fastq)
    'fastqDir',
    # text file with sample names with sample number next to it separated by a space (string, i.e. path/to/libs)    
    'sampleText'
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
    if len(sys.argv)-1 < 3:
       print('need at least 3 arguments: how many cores, where the fastq files are (i.e. working directory), and a text file with the sample names of the fastq files (all on separate line)')
       sys.exit(0)
    cores = sys.argv[1]
    # gather tags from arguments
    wd = sys.argv[2]
    tags = ' ' + inputs[2] + '=' + sys.argv[3]
    # add any optional tags
    tags += addedTags
    # append fastqDir to set working directory
    os.chdir('/rugpfs/fs0/risc_lab/store/npagane/ATACseq') # CHANGE THIS TO FINAL EXECUTABLE DIR
    os.system('snakemake --use-conda --rerun-incomplete --cores ' + cores + ' --directory ' + wd + ' --config' + tags) # CLUSTER CONFIGS HERE
    stop = time.time()
    print('ran took ' + str(1.0*(stop - start)/(60*60)) + ' hours')
