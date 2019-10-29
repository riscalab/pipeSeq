#! usr/bin/env snakemake
## npagane | 190701 | risca lab | fastq2bam pipeline rules

# include this file to incororporate these rules into a Snakefile for execution

import os
import numpy as np
import sys
import helper
import datetime

################################
# parameters and functions
################################

# defaults for parameters set in fastq2bam.py exectuable file

# functions are found in helper.py file

wildcard_constraints:
    post_tag = "\d+"

################################
# commands with custom flags
################################      

align = "(bowtie2 -X2000 -p{threads} -x {config[genomeRef]} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}"

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        r1 = config['fastqDir'] + "/{pre_tag}_R1_{post_tag}.fastq.gz",
        r2 = config['fastqDir'] + "/{pre_tag}_R2_{post_tag}.fastq.gz"
    output:
        config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    params:
        r1 = "{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    run:
        shell("""
              f1=`zcat {input.r1} | awk '$1 ~ /^+/' | wc -l`;
              f2=`zcat {input.r2} | awk '$1 ~ /^+/' | wc -l`;
              l1=`zcat {input.r1} | wc -l`;
              l2=`zcat {input.r2} | wc -l`;
              l14=`expr $l1 '/' 4`;
              l24=`expr $l2 '/' 4`;
              if [ $f1 == $l14 ] && [ $f2 == $l24 ]; then echo 'fast trim'; python {workflow.basedir}/scripts/pyadapter_trimP3V2.py -a {input.r1} -b {input.r2} > {config[sample]}/{wildcards.pre_tag}_adapter_trim.log;
              else echo 'slow trim'; python {workflow.basedir}/scripts/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {config[sample]}/{wildcards.pre_tag}_adapter_trim.log; fi
            """)
        shell("mv {params.r1} {config[sample]}/")
        shell("mv {params.r2} {config[sample]}/")

################################
# QC of fastq files (2)
################################

rule fastqc:
    input:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    output:
        expand(config['sample'] + "/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", "_fastqc.zip", ".fastq"])
    params:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    threads: 2
    run:
        shell("fastqc -t {threads} -o {config[sample]} {input.r1} {input.r2}")
        shell("unpigz {params.r1} {params.r2}")

################################
# align inserts and fastq screen (3)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq",
    output:
        expand(config['sample'] + "/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_screen.html", "_screen.png", "_screen.txt"]),
        bam = config['sample'] + "/{pre_tag}_{post_tag}.trim.unmerged.bam",
        alignLog = config['sample'] + "/{pre_tag}_{post_tag}.trim.alignlog",
    params:
        screen = "screen.log"
    threads: 6
    run:
        shell(align) # align command defined above
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {config[sample]}_{params.screen}")
        shell("mv {config[sample]}_{params.screen} {config[sample]}/{params.screen}")
        shell("mv {wildcards.pre_tag}_R1_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
        shell("mv {wildcards.pre_tag}_R2_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
        shell("pigz {input.unzip1} {input.unzip2}") # zip

################################
# merge sample bams aross lanes (3.5)
################################

rule mergeBamIfNecessary:
    input:
        expand(config['sample'] + "/{{pre_tag}}{lanes}{{post_tag}}.trim.unmerged.bam", lanes=helper.determine_lanes(config["fastqDir"], config["sample"]))
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.bam"
    run:
        if (len(input) == 1):
            shell("echo 'sample {config[sample]} was sequenced in one lane'")
            shell("mv {input} {output}")
        else:
            print('sample {config[sample]} was sequenced in more than one lane; merging BAMs!')
            input_string = ''
            for inp in input:
                input_string += inp + ' '
            shell("samtools merge " + " {output} " + input_string) 

################################
# sort bam files for filtering (4)
################################

rule sortBam:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.bam"
    params:
        "{pre_tag}_{post_tag}.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")

################################
# bam filtering with samtools (5)
################################

rule filterBam:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.all.bam"
    params:
        chrM = config['sample'] + "/{pre_tag}_{post_tag}.trim.st.chrM.bam",
        filterLog = "filtering.log"
    run:
        shell("samtools index {input}")
        shell("echo 'Removing reads from unwanted chromosomes and scaffolds'")
        shell("""
              chrs="";for i in {1..100}; do chrs+=" chr$i"; done; chrs+=" chrX"
              samtools view -b {input} `echo $chrs` > {output}
              """)
        shell("samtools view -b {input} chrM > {params.chrM}")
        shell("echo 'Filtering file {input} by rules filterBam and filter_removeDups_and_enrichTSS' >> {config[sample]}_{params.filterLog}")
        shell("mv {config[sample]}_{params.filterLog} {config[sample]}/{params.filterLog}")

################################
# blacklist and map quality filter, 
# remove duplicates, and TSS enrichment (6)
################################

rule filter_removeDups_and_enrichTSS:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.all.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.bam"
    params:
        filterLog = config['sample'] + "/filtering.log",
        histDupsLog = config['sample'] + "/hist_data_withdups.log",
        histDupsPDF = config['sample'] + "/hist_graphwithdups.pdf",
        dupsLog = config['sample'] + "/dups.log",
        histNoDupsLog = config['sample'] + "/hist_data_withoutdups.log",
        histNoDupsPDF = config['sample'] + "/hist_graphwithoutdups.pdf"
    run:
        # filter by blacklist if provided
        if os.path.exists(config["blacklist"]):
            shell("echo 'Removing blacklisted reads'")
            shell("bedtools intersect -v -abam {input} -b {config[blacklist]} -wa > {config[sample]}_temp.bam") # produces temp file
            ftp = "{config[sample]}/{wildcards.pre_tag}_{wildcards.post_tag}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 {config[sample]}_temp.bam -o " + ftp)
            shell("rm {config[sample]}_temp.bam") # remove temp file
            shell("echo 'Blacklist filtered using file {config[blacklist]}.' >> {params.filterLog}")
        else:
            ftp = input
            shell("echo 'Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")
        # filter by map quality if provided
        if int(config["mapq"]) > 0:
            shell("echo 'Removing low quality reads'")
            tmp = str(ftp).split('.bam')[0] + '.qft.bam'
            shell("samtools view -bh -f 0x2 -q {config[mapq]} " + str(ftp) + " -o " + tmp)
            ftp = tmp
            shell("echo 'Filtered with mapping quality filter {config[mapq]}.' >> {params.filterLog}")
        else:
            shell("echo 'Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")
        # histogram with duplicates
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=1000 STOP_AFTER=5000000")
        # remove duplicates
        shell("echo 'Removing duplicates'")
        tmp = ftp.split('.bam')[0] + '.rmdup.bam'
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histNoDupsLog} H={params.histNoDupsPDF} W=1000 STOP_AFTER=5000000")
        # TSS enrichment if provided
        if os.path.exists(config["TSS"]):
            tmp = ftp + '.RefSeqTSS'
            shell(workflow.basedir + "/scripts/pyMakeVplot_css_v01.py -a " + ftp + " -b {config[TSS]} -e 2000 -p ends -s 5 -v -u --atac -o " + tmp)
        else:
            shell("echo 'TSS BED file not provided. not creating TSS enrichment profile'")
        # cleanup directory
        shell("if [ ! -d {config[sample]}/00_source ]; then mkdir {config[sample]}/00_source; fi") 
        shell("mv {config[sample]}/*.all.bam {config[sample]}/00_source/")
        shell("mv {config[sample]}/*.chrM.bam {config[sample]}/00_source/")

