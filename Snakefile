#! usr/bin/env snakemake
# npagane | 190701 | risca lab | snakefile including main rules and target rule

import .rules.helper as helper

################################
# pipeline output check
################################

rule all:
    input:
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"
            ), config['fastqDir'], config['sample']
        )

################################
# rules for fastq2bam
################################

include: "rules/fastq2bam.smk"

