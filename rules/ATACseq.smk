# usr/bin/env snakemake
## npagane | 190701 | risca lab | ATACseq pipeline rules

# include this file to incororporate these rules into a Snakefile for execution

import os
import sys
import datetime
import helper

################################
# parameters and functions
################################

# defaults for parameters set in ATACseq.py exectuable file

# functions are found in helper.py file

wildcard_constraints:
    post_tag = "\d+"

################################
# commands with custom flags
################################      

callpeak = helper.callpeak # go to helper file to see / edit callpeak command
bdgcmpfc = "macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m FE" # fold enrichment 
preS = """calc(){{ awk "BEGIN {{ print "$*" }}"; }}; num=`wc -l {input.bed} | awk '{{print $1}}'`; den=1000000;"""
S = """S=`calc $num/$den | awk '{{printf("%i", $1)}}'`;""" 
bdgcmppval = "macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m ppois -S $S" # p value
bam2bg = helper.bam2bg # go to helper file to see / edit bam2bg command

################################
# align at insertion center (1)
################################

rule ATACoffset:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.bam"
    params:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.temp.bam"
    threads: 6
    run:
        shell("samtools index {input}") # suppress the pysam/htslib warning about the index file
        shell("alignmentSieve --numberOfProcessors {threads} --ATACshift --bam {input} -o {params}")
        shell("samtools sort -O bam -o {output} {params}") #sort (dont use picard it is too strict about bam formatting)
        shell("samtools index {output}") # regenerate index file
        shell("rm {params}")

################################
# call peaks to make bed (2)
################################

rule callPeakSummits:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.atac.bam"
    output:
        expand(config['sample'] + "/peakCalls/{{pre_tag}}_{{post_tag}}{{ext,.*}}.rmdup.atac{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak", "_treat_pileup.bdg", "_control_lambda.bdg"])
    params:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac"
    conda:
        "../envs/macs2_python2.yml" # path relative to current file, not working directory
    shell:
        callpeak # macs2 callpeak command defined above

################################
# make tracks genomecov (3a)
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

################################
# signal generation to form bedgraph
# with fold change values (3b)
################################

rule makeBedGraphSignalFC:
    input:
        treat = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_treat_pileup.bdg",
        control = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_control_lambda.bdg"
    output:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac_FE.bdg"
    conda:
        "../envs/macs2_python2.yml" # path relative to current file, not working directory
    shell:
        bdgcmpfc # macs2 bgdcmp command defined above for fold change

################################
# signal generation to form bedgraph
# with fold change values (4b)
################################

rule makeBedGraphSignalPval:
    input:
        treat = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_treat_pileup.bdg",
        control = config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext}.rmdup.atac_control_lambda.bdg",
        bed = config['sample'] + "/tracks/{pre_tag}_{post_tag}{ext}.rmdup.atac.bed"
    output:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac_pval.bdg"
    conda:
        "../envs/macs2_python2.yml" # path relative to current file, not working directory
    shell:
        preS + ' ' + S + ' ' + bdgcmppval # macs2 bgdcmp command defined above for p value

################################
# create bigwig for tracks (5)
################################

rule bedGraph2bigWig:
    input:
        config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext}.rmdup.atac{ext2}.bdg"
    output:
        bw = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bw"
    params:
        qft = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bdg.clip",
        st = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bdg.st.clip"
    run:
        shell("bedtools slop -i {input} -g {config[chromSize]} -b 0 | bedClip stdin {config[chromSize]} {params.qft}")
        shell("LC_COLLATE=C sort -k1,1 -k2,2n {params.qft} > {params.st}")
        shell("bedGraphToBigWig {params.st} {config[chromSize]} {output.bw}")
        shell("rm {params.qft} {params.st}")
        # cleanup directory
        shell("if [ ! -d {config[sample]}/00_source ]; then mkdir {config[sample]}/00_source; fi")
        shell("mv {config[sample]}/*atac.bam {config[sample]}/00_source/")
        shell("mv {config[sample]}/*bai {config[sample]}/00_source/")
