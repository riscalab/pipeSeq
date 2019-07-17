#! usr/bin/env snakemake
# npagane | 190701 | risca lab | snakefile including main rules and target rule

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
