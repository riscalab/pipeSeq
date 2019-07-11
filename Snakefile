#! /bin/bash snakemake
# npagane | 190618 | risca lab | snakefile for fastq2bam pipeline

################################
# pipeline output check
################################

rule all:
    input:
        "fastq2bamRunSummary.log"

################################
# rules for fastq2bam
################################

include: "./rules/fastq2bam.smk"

