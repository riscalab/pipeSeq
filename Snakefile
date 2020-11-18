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
            # all tracks and peaks
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
            ),
            # and TSS enrichment files
            helper.customFileExpand(
                helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                    ".trim.st.all.blft.qft.RefSeqTSS.log",
                    ".trim.st.all.qft.RefSeqTSS.log",
                    ".trim.st.all.blft.RefSeqTSS.log",
                    ".trim.st.all.RefSeqTSS.log"
                ), config['fastqDir'], config['sample'], 
            )
elif config['pipe'] == 'CUTnTag':
    # THIS WHILL CHANGE IN THE FUTURE. WILL ADD PEAK STUFF BUT FOR NOW MAKE IT LIKE fastq2bam
    rule all:
        input:
            helper.customFileExpand(
                # final BAMs
                helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                    ".trim.st.all.blft.qft.rmdup.bam",
                    ".trim.st.all.qft.rmdup.bam",
                    ".trim.st.all.blft.rmdup.bam",
                    ".trim.st.all.rmdup.bam"
                ), config['fastqDir'], config['sample']
            )
 
else: # default to fastq2bam
    rule all:
        input:
            helper.customFileExpand(
                # final BAMs
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

if config['pipe'] == 'fastq2bam':
    include: "rules/fastq2bam/trimAdapters.smk"
    include: "rules/fastq2bam/fastqcAndScreen.smk"
    include: "rules/fastq2bam/alignInserts.smk"
    include: "rules/fastq2bam/mergeBamIfNecessary.smk"
    include: "rules/fastq2bam/sortBam.smk"
    include: "rules/fastq2bam/filterBam.smk"
    include: "rules/fastq2bam/blacklistFilter.smk"
    include: "rules/fastq2bam/qualityFilter.smk"
    include: "rules/fastq2bam/removeDuplicates.smk"

if config['pipe'] == 'ATACseq':
    # fastq2bam stuff
    include: "rules/fastq2bam/trimAdapters.smk"
    include: "rules/fastq2bam/fastqcAndScreen.smk"
    include: "rules/fastq2bam/alignInserts.smk"
    include: "rules/fastq2bam/mergeBamIfNecessary.smk"
    include: "rules/fastq2bam/sortBam.smk"
    include: "rules/fastq2bam/filterBam.smk"
    include: "rules/fastq2bam/blacklistFilter.smk"
    include: "rules/fastq2bam/qualityFilter.smk"
    include: "rules/fastq2bam/removeDuplicates.smk"
    # ATACseq stuff
    include: "rules/ATACseq/ATACoffset.smk"
    include: "rules/ATACseq/enrichTSS.smk"
    include: "rules/ATACseq/callPeakSummits.smk"
    include: "rules/ATACseq/bam2bed.smk"
    include: "rules/ATACseq/makeBedGraphSignalFC.smk"
    include: "rules/ATACseq/makeBedGraphSignalPval.smk"
    include: "rules/ATACseq/bedGraph2bigWig.smk"

if config['pipe'] == 'CUTnTag':
    # fastq2bam stuff
    include: "rules/fastq2bam/trimAdapters.smk"
    include: "rules/fastq2bam/fastqcAndScreen.smk"
    include: "rules/CUTnTag/alignInserts.smk" # only aligning is different for now
    include: "rules/fastq2bam/mergeBamIfNecessary.smk"
    include: "rules/fastq2bam/sortBam.smk"
    include: "rules/fastq2bam/filterBam.smk"
    include: "rules/fastq2bam/blacklistFilter.smk"
    include: "rules/fastq2bam/qualityFilter.smk"
    include: "rules/fastq2bam/removeDuplicates.smk"

