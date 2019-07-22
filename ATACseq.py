#! usr/bin/env python
# npagane | 190701 | risca lab | execute snakemake for ATACseq

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
    # genome reference for alignment (string, i.e. path/to/genome)
    '--genomeRef': "/rugpfs/fs0/risc_lab/scratch/nvelez/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome",
    # blacklist for filtering (string, i.e. path/to/blacklist)
    '--blacklist': "/rugpfs/fs0/risc_lab/scratch/nvelez/blacklists/ATAC_blacklist.bed",
    # TSS for that genome (string, i.e. path/to/TSS)
    '--TSS': "/rugpfs/fs0/risc_lab/store/vrisca/lab-shared/dl-annotations/hg38/GENCODE/gencode.v30.basic.annotation.tss.bed",
    # the map quality threshold for alignment (int, i.e. 30)
    '--mapq': '30', 
    # whether there are fastq files for the index reads or not (string, i.e. True)
    '--index': 'False',
    # exclusion file to exclude certain samples from analysis (string, i.e. path/to/exclude)
    '--exclude': '',
    # any snakemake flags for compilation (string, i.e. unlock)
    '--snakemake': '', 
    # sizes of chromosomes (string, i.e. path/to/genome/file)
    '--chromSize': "/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes",
   
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
    snakeFlags = ''
    for tag in optTags.keys():
        tempTag = ' '.join(sys.argv).split(tag)
        if len(tempTag) == 1:
            addedTags += " '" + tag[2:] + '="' + optTags[tag] + '"' + "'"
        else:
            tempFlag = tempTag[1].split(' ')[1]
            if tag == "--snakemake":
                snakeFlags += " --" + tempFlag
            else:
                addedTags += " '" + tag[2:] + '="' + tempFlag + '"' + "'"
    # set working directory to project directory and figure out where the snakefile is
    os.chdir(wd)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    # execute snakemake and pass path to Snakefile with proper configs
    os.system('snakemake --snakefile ' + dir_path +  '/Snakefile --use-conda --conda-prefix ./.snakemake --rerun-incomplete --cores ' + cores + ' ' + snakeFlags + ' --config' + addedTags) 
    stop = time.time()
    print('ran took ' + str(1.0*(stop - start)/(60*60)) + ' hours')
