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
                ".trim.st.all.blft.qft.rmdup.atac.bw",
                ".trim.st.all.qft.rmdup.atac.bw",
                ".trim.st.all.blft.rmdup.atac.bw",
                ".trim.st.all.rmdup.atac.bw"
            ), config['fastqDir'], config['sample'], 'tracks'
        ),
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.atac_FE.bw",
                ".trim.st.all.qft.rmdup.atac_FE.bw",
                ".trim.st.all.blft.rmdup.atac_FE.bw",
                ".trim.st.all.rmdup.atac_FE.bw"
            ), config['fastqDir'], config['sample'], 'peakCalls'
        ),
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.atac_pval.bw",
                ".trim.st.all.qft.rmdup.atac_pval.bw",
                ".trim.st.all.blft.rmdup.atac_pval.bw",
                ".trim.st.all.rmdup.atac_pval.bw"
            ), config['fastqDir'], config['sample'], 'peakCalls'
        )

################################
# rules for ATACseq
################################

include: "rules/ATACseq.smk"
