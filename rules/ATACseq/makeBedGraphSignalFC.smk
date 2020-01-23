#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

################################
# commands with custom flags
################################

bdgcmpfc = "macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m FE" # fold enrichment 

################################
# signal generation to form bedgraph
# with fold change values (4a)
################################

rule makeBedGraphSignalFC:
    input:
        treat = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_treat_pileup.bdg",
        control = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_control_lambda.bdg"
    output:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac_FE.bdg"
    conda:
        "../../envs/macs2_python2.yaml" # path relative to current file, not working directory
    shell:
        bdgcmpfc # macs2 bgdcmp command defined above for fold change
