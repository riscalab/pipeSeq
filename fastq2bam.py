#! /bin/bash python3
# npagane | 190618 | risca lab | execute fastq2bam pipeline remotely

import os
import sys


# parameters needed to insert as argruments to execute
inputs = [
    'core',
    'fastqDir',
    'sampleText',
    'genomeRef',
    'blacklist',
    'fastaRef'
]

if __name__ == '__main__':
    if len(sys.argv)-1 < 3:
       print('need at least 3 arguments: how many cores, where the fastq files are (i.e. working directory), and a text file with the sample names of the fastq files (all on separate line)')
       sys.exit(0)
    cores = sys.argv[1].strip()
    # gather tags from arguments
    tags = ''
    wd = ' ' + sys.argv[2].strip()
    for i in range(2, len(inputs)):
        if i < (len(sys.argv) - 1):
            tags += ' ' + str(inputs[i]) + '=' + sys.argv[i+1].strip()
        else:
            tags += ' ' + str(inputs[i]) + "=''"
    # append fastqDir to set working directory
    os.chdir('/rugpfs/fs0/risc_lab/store/npagane/fastq2bam') # CHANGE THIS TO FINAL EXECUTABLE DIR
    os.system('snakemake --cores ' + cores + ' --directory' + wd + ' --config' + tags) # CLUSTER CONFIGS HERE
    # add -n tag for dry run # CHANGE CORE NUMBER
