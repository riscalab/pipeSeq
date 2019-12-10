#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

import sys
sys.path.append('{workflow.basedir}')
import rules.helper

################################
# commands with custom flags
################################      

bam2bg = helper.bam2bg # go to helper file to see / edit bam2bg command

################################
# make tracks genomecov (3)
################################

rule bam2bed:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.atac.bam"
    output:
        bg = config['sample'] + "/tracks/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.bdg",
        bed = config['sample'] + "/tracks/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.bed"
    run:
        shell("bedtools bamtobed -i {input} > {output.bed}") # bam to bed necessary for bigwig analysis
        shell(bam2bg) # bam to bedgraph command defined above
