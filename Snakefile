#! /bin/bash snakemake

################################
# pipeline output check
################################

rule all:
    input:
        "ATACseqRunSummary.log"

################################
# rules for fastq2bam
################################

include: "rules/fastq2bam.smk"

################################
# rules for ATACseq
################################

include: "rules/ATACseq.smk"
