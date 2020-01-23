#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

################################
# commands with custom flags
################################  

preS = """calc(){{ awk "BEGIN {{ print "$*" }}"; }}; num=`wc -l {input.bed} | awk '{{print $1}}'`; den=1000000;"""
S = """S=`calc $num/$den | awk '{{printf("%i", $1)}}'`;"""
bdgcmppval = "macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m ppois -S $S" # p value

################################
# signal generation to form bedgraph
# from significantly enriched values (4b)
################################

rule makeBedGraphSignalPval:
    input:
        treat = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_treat_pileup.bdg",
        control = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_control_lambda.bdg",
        bed = config['sample'] + "/tracks/{pre_tag}_{post_tag}{ext}.rmdup.atac.bed"
    output:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac_pval.bdg"
    conda:
        "../../envs/macs2_python2.yaml" # path relative to current file, not working directory
    shell:
        preS + ' ' + S + ' ' + bdgcmppval # macs2 bgdcmp command defined above for p value
