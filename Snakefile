#! /bin/bash snakemake

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

