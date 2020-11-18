#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

import sys
sys.path.append('{workflow.basedir}')
import rules.helper

################################
# commands with custom flags
################################      

align = helper.align_CUTnTag # go to helper file to see / edit align command

################################
# align inserts to genome (3)
################################

rule alignInserts:
    input:
        unzip1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq",
    output:
        bam = config['sample'] + "/{pre_tag}_{post_tag}.trim.unmerged.bam",
        alignLog = config['sample'] + "/{pre_tag}_{post_tag}.trim.alignlog",
    threads: 9
    run:
        shell(align) # align command defined above
        shell("pigz {input.unzip1} {input.unzip2}") # zip
