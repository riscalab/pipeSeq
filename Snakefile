#! /bin/bash snakemake
# npagane | 190701 | risca lab | snakefile for ATACseq pipeline

import os
import sys
import datetime

################################
# pipeline output check
################################

rule all:
    input:
        "ATACseqRunSummary.log"

################################
# rules for fastq2bam
################################

include: "./../fastq2bam/rules/fastq2bam.smk"

################################
# rules for ATACseq
################################

include: "./rules/ATACseq.smk"
