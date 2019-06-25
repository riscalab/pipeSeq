#! /bin/bash python3
# npagane | 190618 | risca lab | snakefile for fastq2bam pipeline

################################
# parameters and functions
################################

import os

# reassign the file to sample names to a list of the actual samples
SAMPLES = []
with open(config['sampleText']) as ftp:
    lines = ftp.readlines()
    lane = {}
    cnt = 1
    for line in lines:
        sample = line.strip()
        SAMPLES.append(sample)
        lane[sample] = str(cnt)
        cnt += 1

# additional fastq file name variables like runs etc. 
TAGS = ['_R1', '_R2', '_I1', '_I2'] 

# customized expand function to be compatible with a python dictionary
def customSeqFileExpand(iden, ext, wd = False): # enter wildcards argument if it starts working
    strout = []
    for sample in SAMPLES:
        # see if just the name of the files or with the sample directories
        if wd:
            primer = sample + '/'
        else:
            primer = ''
        # create files
        for i in iden:
            ftp = primer + sample + '_S' + lane[sample] + i + '_001' + ext 
            strout.append(ftp)
    return strout

# customized expand function to be compatible with a python dictionary
def customFileExpand(out): # enter wildcards argument if it starts working
    strout = []
    for sample in SAMPLES:
        ftp = sample + '/' + out
        strout.append(ftp)
    return strout

# conditional expand function upon one condition for inputs/output
def conditionalExpand_1(condition, true, false):
    if condition:
        fpt = true
    else:
        fpt = false
    return fpt

# conditional expand function upon 2 conditions for inputs/output
def conditionalExpand_2(condition1, condition2, truetrue, truefalse, falsetrue, falsefalse):
    if condition1:
        if condition2:
            fpt = truetrue
        else:
            fpt = truefalse
    else:
        if condition2:
            fpt = falsetrue
        else:
            fpt = falsetrue
    return fpt

# be careful of set variable (usually 001); somehow make wildcard?

# map quality
MAPQ = 30 # maybe make this a prompt at somepoint.. but for now, this is fine

################################
# pipeline output check
################################

rule all:
    input:
        customSeqFileExpand(TAGS, '.fastq.gz', True), # path rule 1
        customSeqFileExpand(['_R1', '_R2'], '_fastqc.html', True), # path rule 2
        customSeqFileExpand(['_R1', '_R2'], '.trim.fastq.gz', True), # path rule 2
        customSeqFileExpand([''], '.trim.bam', True), # path rule 3
        customSeqFileExpand([''], '.trim.st.bam', True), # path rule 4
        customSeqFileExpand([''], '.trim.st.all.bam', True), # path rule 5
         customFileExpand("/hist_data_withdups.log"), #path rule 8
        customFileExpand("hist_graphwithoutdups.pdf") # path rule 9

################################
# organize sample directories (1)
################################

rule createSampleDirectories:
    input:
        "{sample}_{lane}_{tag}_001.fastq.gz"
    output:
        "{sample}/{sample}_{lane}_{tag}_001.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# QC of fastq files and trim (2)
################################

rule fastqQC_and_trimAdapters:
    input:
        r1 = "{sample}/{sample}_{lane}_R1_001.fastq.gz",
        r2 = "{sample}/{sample}_{lane}_R2_001.fastq.gz"
    output:
        "{sample}/{sample}_{lane}_R1_001_fastqc.html",
        "{sample}/{sample}_{lane}_R2_001_fastqc.html",
        #"{sample}/{sample}_{lane}_R1_001_fastqc.zip", zip files for html
        #"{sample}/{sample}_{lane}_R2_001_fastqc.zip", zip file for html
        "{sample}/{sample}_{lane}_R1_001.trim.fastq.gz",
        "{sample}/{sample}_{lane}_R2_001.trim.fastq.gz"
    run:
        shell("fastqc -o {wildcards.sample} {input.r1} {input.r2}")
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")

################################
# align inserts and fastq screen (3)
################################

rule alignInserts_and_fastqScreen:
    input:
        zip1 = "{sample}/{sample}_{lane}_R1_001.trim.fastq.gz",
        zip2 = "{sample}/{sample}_{lane}_R2_001.trim.fastq.gz",
    output:
        bam = "{sample}/{sample}_{lane}_001.trim.bam",
        alignLog = "{sample}/{sample}_{lane}_001.trim.alignlog",
    params:
        unzip1 = "{sample}/{sample}_{lane}_R1_001.trim.fastq",
        unzip2 = "{sample}/{sample}_{lane}_R2_001.trim.fastq",
        ref = config['genomeRef'],
        screen = "screen.log"
    run:
        shell("gunzip {input.zip1}") # unzip
        shell("gunzip {input.zip2}") # unzip
        shell("(bowtie2 -p28 -x {params.ref} -1 {params.unzip1} -2 {params.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}")
        shell("fastq_screen --aligner bowtie2 {params.unzip1} {params.unzip2}")
        shell("{params.screen}")
        shell("gzip {output.unzip1}") # zip
        shell("gzip {output.unzip2}") # zip

################################
# sort bam files for filtering (4)
################################

rule sortBam:
    input:
        "{sample}/{sample}_{lane}_001.trim.bam"
    output:
        "{sample}/{sample}_{lane}_001.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={ouput}  SORT_ORDER=coordinate")


################################
# bam filtering with samtools (5)
################################

rule filterBam:
    input:
        "{sample}/{sample}_{lane}_001.trim.st.bam"
    output:
        allBam = "{sample}/{sample}_{lane}_001.trim.st.all.bam",
        filterLog = "{sample}/filtering.log"
    params:
        chrs = "`samtools view -H *.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v _gl | grep -v Y | grep -v hap | grep -v random | grep -v v1 | grep -v v2`"
    run:
        shell("samtools index {input}")
        shell("echo 'sh02a_filter_bam.sh: Removing reads from unwanted chromosomes and scaffolds'")
        shell("samtools view -b {input} `echo {params.chrs}` > {output.allBam}")
        shell("echo 'Filtering file {input} by script sh02a_filter_bam.sh' >> {output.filterLog}")


################################
# filter by blacklists (6)
################################

rule filterBamWithBlacklist:
    input:
        "{sample}/{sample}_{lane}_001.trim.st.all.bam"
    output:
        bfltBam = conditionalExpand_1(os.path.exists(config['blacklist']), 
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam", 
            "{sample}/{sample}_{lane}_001.trim.st.all.bam"),
        filterLog = "{sample}/filtering.log"
    params:
        blacklist = config['blacklist']
    run:
        # filter by blacklist if provided
        if os.path.exists(params.blacklist):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > temp.bam") # produces temp file
            shell("samtools view -bh -f 0x2 temp.bam -o {output.blftBam}")
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {output.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {output.filterLog}")

################################
# filter by map quality (7)
################################

rule filterBamWithMapQ:
    input:
        inFile = conditionalExpand_1(MAPQ, 
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam", 
            "{sample}/{sample}_{lane}_001.trim.st.all.bam")
    output:
        outFile = conditionalExpand_1(MAPQ, 
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.bam", 
            "{sample}/{sample}_{lane}_001.trim.st.all.qft.bam"),
        filterLog = "{sample}/filtering.log"
    params:
        mapq = MAPQ
    run:
        if MAPQ > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            shell("samtools view -bh -f 0x2 -q {params.mapq} {input.inFile} -o {output.outFile}")
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {output.filterLog}")
            shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {output.filterLog}")

################################
# histogram with duplicates (8)
################################

rule histogramWithDuplicates_and_removeDuplicates:
    input:
        inFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.qft.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.bam")
    output:
        outFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.qft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.rmdup.bam"),
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log"
    params:
        blacklist = config['blacklist'],
        mapq = MAPQ
    run:
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I={input.inFile} O={output.histDupsLog} H={output.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        shell("picard MarkDuplicates I={input.inFile} O={output.outFile} METRICS_FILE={output.dupsLog} REMOVE_DUPLICATES=true")
        shell("samtools index {output.outFile}")

################################
# histogram with removed dups (9)
################################

rule histogramWithNoDuplicates:
    input:
        inFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.qft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.blft.rmdup.bam",
            "{sample}/{sample}_{lane}_001.trim.st.all.rmdup.bam")
    output:
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    params:
        blacklist = config['blacklist'],
        mapq = 30
    run:
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I={input.inFile} O={output.histNoDupsLog} H={output.histNoDupsPDF} W=1000 STOP_AFTER=50000")


################################
# clean up ??
################################
'''
rule cleanUp:
    input:
         "{sample}/hist_data_withoutdups.log"
    output:
        "{sample}/00_source/"

'''