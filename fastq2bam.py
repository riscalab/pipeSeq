#! /bin/bash python3
# npagane | 190618 | risca lab | execute fastq2bam pipeline remotely

import os

# initialize dictionary to store parameters
keys = ['sampleText', 'fastqDir', 'genomeRef', 'blacklist', 'fastaRef']

prompts = {
    'fastqDir': 'global path to directory with fastq.gz files (all commands will be executed here): ',
    'sampleText': 'global path to text file with sample names: ',
    'genomeRef': 'global path to genome index  reference file (optional): ',
    'blacklist': 'global path to the genome reference for chr sizes (optional): ',
    'fastaRef' : 'global path to the FASTA genome reference file (optional): '
}

if __name__ == '__main__':
    # gather tags from standard output
    tags = ''
    wd = ' ' + input(prompts["fastqDir"]).strip()
    for i in keys:
        if i != "fastqDir":
            tags += ' ' + str(i) + '=' + input(prompts[i]).strip()
    # append fastqDir to set working directory
    os.chdir('/rugpfs/fs0/risc_lab/store/npagane/fastq2bam') # CHANGE THIS TO FINAL EXECUTABLE DIR
    os.system('snakemake --directory' + wd + ' --config' + tags) # CLUSTER CONFIGS HERE
    # add -n tag for dry run