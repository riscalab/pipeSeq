#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

import sys
sys.path.append('../')
import helper

################################
# commands with custom flags
################################      

align = helper.align # go to helper file to see / edit align command

################################
# align inserts and fastq screen (3)
# ##############################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq",
    output:
        expand(config['sample'] + "/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_screen.html", "_screen.png", "_screen.txt"]),
        bam = config['sample'] + "/{pre_tag}_{post_tag}.trim.unmerged.bam",
        alignLog = config['sample'] + "/{pre_tag}_{post_tag}.trim.alignlog",
    threads: 18
    run:
        shell(align) # align command defined above
        # go to /rugpfs/fs0/risc_lab/store/risc_soft/anaconda2/envs/fastq2bam/share/fastq-screen-0.13.0-0/fastq_screen.conf to make changes to fastq_screen
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}")
        shell("mv {wildcards.pre_tag}_R1_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
        shell("mv {wildcards.pre_tag}_R2_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
        shell("pigz {input.unzip1} {input.unzip2}") # zip
