#! usr/bin/env snakemake
## npagane | risca lab | snakefile including main rules and target rule

import rules.helper as helper

wildcard_constraints:
    post_tag = "\d+"

################################
# pipeline output check
################################

if config['pipe'] == 'ATACseq':
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
else: # default to fastq2bam
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
# rules for pipeline
################################

# assume general fastq2bam rules (THIS WILL NEED TO CHANGE WHEN RICC SEQ ADDED)

include: "rules/fastq2bam/trimAdapters.smk"
include: "rules/fastq2bam/fastqc.smk"
include: "rules/fastq2bam/alignInserts_and_fastqScreen.smk"
include: "rules/fastq2bam/mergeBamIfNecessary.smk"
include: "rules/fastq2bam/sortBam.smk"
include: "rules/fastq2bam/filterBam.smk"
include: "rules/fastq2bam/blacklistFilter_removeDups_and_enrichTSS.smk"

if config['pipe'] == 'ATACseq':
    include: "rules/ATACseq/ATACoffset.smk"
    include: "rules/ATACseq/callPeakSummits.smk"
    include: "rules/ATACseq/bam2bed.smk"
    include: "rules/ATACseq/makeBedGraphSignalFC.smk"
    include: "rules/ATACseq/makeBedGraphSignalPval.smk"
    include: "rules/ATACseq/bedGraph2bigWig.smk"
