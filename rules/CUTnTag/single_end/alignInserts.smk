#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

import sys
sys.path.append('{workflow.basedir}')
import rules.helper

################################
# commands with custom flags
################################      

align = helper.align_CUTnTag_SE # go to helper file to see / edit align command

################################
# align inserts to genome (3)
################################

rule alignInserts:
    input:
        unzip = config['sample'] + "/{pre_tag}_{post_tag}.fastq",
    output:
        bam = config['sample'] + "/{pre_tag}_{post_tag}.unmerged.bam",
        alignLog = config['sample'] + "/{pre_tag}_{post_tag}.alignlog",
    threads: 9
    run:
        shell(align) # align command defined above
        shell("pigz {input.unzip}") # zip
