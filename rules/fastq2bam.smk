#! fastq2bam pipeline rules

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

# determine if there are index fastq files
if config['index'] == 'True':
    TAGS = ['R1', 'R2', 'I1', 'I2'] 
else:
    TAGS = ['R1', 'R2']

wildcard_constraints:
    post_tag = "\d+"

################################
# organize sample directories (0)
################################

rule createSampleDirectories:
    input:
        "{pre_tag}_{tag}_{post_tag}.fastq.gz"
    output:
        "{sample}/{pre_tag}_{tag}_{post_tag}.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        expand("{{sample}}/{{pre_tag}}_{tag}_{{post_tag}}.fastq.gz", tag=TAGS), # move all files to sample dirs
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.fastq.gz"
    output:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    params:
        r1 = "{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    run:
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")
        shell("mv {params.r1} {wildcards.sample}/")
        shell("mv {params.r2} {wildcards.sample}/")

################################
# QC of fastq files (2)
################################

rule fastqQC:
    input:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    output:
        expand("{{sample}}/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", "_fastqc.zip", ".fastq"])
    params:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    run:
        shell("fastqc -o {wildcards.sample} {input.r1} {input.r2}")
        shell("unpigz {params.r1} {params.r2}")

################################
# align inserts and fastq screen (3)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq",
    output:
        expand("{{sample}}/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_screen.html", "_screen.png", "_screen.txt"]),
        bam = "{sample}/{pre_tag}_{post_tag}.trim.bam",
        alignLog = "{sample}/{pre_tag}_{post_tag}.trim.alignlog"
    params:
        ref = config['genomeRef'],
        screen = "screen.log"
    run:
        shell("(bowtie2 -p18 -x {params.ref} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}")
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {wildcards.sample}_{params.screen}")
        shell("mv {wildcards.sample}_{params.screen} {wildcards.sample}/{params.screen}")
        shell("mv {wildcards.pre_tag}_R1_{wildcards.post_tag}.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.pre_tag}_R2_{wildcards.post_tag}.trim_screen.* {wildcards.sample}/")
        shell("pigz {input.unzip1} {input.unzip2}") # zip

################################
# sort bam files for filtering (4)
################################

rule sortBam:
    input:
        "{sample}/{pre_tag}_{post_tag}.trim.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.trim.st.bam"
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
        "{sample}/{pre_tag}_{post_tag}.trim.st.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.trim.st.all.bam"
    params:
        chrs = "samtools view -H {sample}/{pre_tag}_{post_tag}.trim.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v _gl | grep -v Y | grep -v hap | grep -v random | grep -v v1 | grep -v v2",
        filterLog = "filtering.log"
    run:
        shell("samtools index {input}")
        shell("echo 'sh02a_filter_bam.sh: Removing reads from unwanted chromosomes and scaffolds'")
        shell("samtools view -b {input} `echo {params.chrs}` > {output}")
        shell("echo 'Filtering file {input} by script sh02a_filter_bam.sh' >> {wildcards.sample}_{params.filterLog}")
        shell("mv {wildcards.sample}_{params.filterLog} {wildcards.sample}/{params.filterLog}")

################################
# filter and remove duplicates (6)
################################

rule filter_and_removeDuplicates:
    input:
        "{sample}/{pre_tag}_{post_tag}.trim.st.all.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.bam"
    params:
        blacklist = config['blacklist'],
        mapq = int(config['mapq']),
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log",
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    run:
        # filter by blacklist if provided
        if os.path.exists(params.blacklist):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads'")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > {wildcards.sample}_temp.bam") # produces temp file
            ftp = "{wildcards.sample}/{wildcards.pre_tag}_{wildcards.post_tag}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 {wildcards.sample}_temp.bam -o " + ftp)
            shell("rm {wildcards.sample}_temp.bam") # remove temp file
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {params.filterLog}")
        else:
            ftp = input
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")
        # filter by map quality if provided
        if params.mapq > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            tmp = str(ftp).split('.bam')[0] + '.qft.bam'
            shell("samtools view -bh -f 0x2 -q {params.mapq} " + str(ftp) + " -o " + tmp)
            ftp = tmp
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")
        # histogram with duplicates
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=600 STOP_AFTER=5000000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        tmp = ftp.split('.bam')[0] + '.rmdup.bam'
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histNoDupsLog} H={params.histNoDupsPDF} W=600 STOP_AFTER=5000000")
        # cleanup directory
        shell("mkdir {wildcards.sample}/00_source")
        shell("mv {wildcards.sample}/*.all.bam {wildcards.sample}/00_source/")
        shell("mv {wildcards.sample}/*.st.bam.bai {wildcards.sample}/00_source/")

################################
# success and summary (6)
################################

rule fastq2bamSummary:
    input:
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"
            ), config['exclude']
        )
    output:
        "fastq2bamRunSummary.log"
    run:
        # make nice pdf of insert distributions
        FILES = np.unique(list(helper.findFiles(config['exclude']).keys())) # order the samples
        libs = "libs.txt"
        with open(libs, "w") as f:
            for i in FILES:
                f.write(i + '\n')
        shell("Rscript " + workflow.basedir + "/scripts/plotisds_v2.R " + libs + " hist_data_withoutdups")
        shell("rm " + libs)
        print('\n###########################')
        print('fastq2bam pipeline complete')
        print('\n###########################')
        with open(output[0], "w") as f:
            f.write('user: ' + os.environ.get('USER') + '\n')
            f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
            f.write('env: fastq2bam\n')
            f.write("SOFTWARE\n")
            f.write("########\n")
            f.write("python version: " + str(sys.version_info[0]) + '\n')
            f.write("pyadapter_trim version: python3 compatible (v1)" + '\n')
            f.write("fastqc version: " + os.popen("fastqc --version").read() + '\n')
            f.write("bowtie2 version: " + os.popen("bowtie2 --version").read() + '\n')
            f.write("samtools version: " + os.popen("samtools --version").read() + '\n')
            f.write("picard version: 2.20.2-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
            f.write("bedtools version: " + os.popen("bedtools --version").read() + '\n')
            f.write("\n\n")
            f.write("PARAMETERS\n")
            f.write("##########\n")
            f.write("genome reference for aligning: " + config["genomeRef"] + '\n')
            f.write("blacklist for filtering: " + config["blacklist"] + '\n')
            f.write("map quality threshold for filtering: " + config["mapq"] + '\n')
            f.write("align command: (bowtie2 -p28 -x {genomeReference} -1 {input.R1} -2 {input.R2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}" + '\n')
            f.write("SUMMARY\n")
            f.write("#######\n")
