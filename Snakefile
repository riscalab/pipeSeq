#! /bin/bash python3
# npagane | 190618 | risca lab | snakefile for fastq2bam pipeline

################################
# parameters and configurations
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
def customExpand(iden, ext, wd = False): # enter wildcards argument if it starts working
    strout = []
    for sample in SAMPLES:
        # see if just the name of the files or with the sample directories
        if wd:
            primer = sample + '/'
        else:
            primer = ''
        for i in iden:
            ftp = primer + sample + '_S' + lane[sample] + i + '_001' + ext
            strout.append(ftp)
    return strout

# be careful of set variable (usually 001); somehow make wildcard?

# map quality??
MAPQ = 30

################################
# pipeline output check
################################

rule all:
    input:
        customExpand(TAGS, '.fastq.gz', True), # path rule 1
        customExpand(['_R1', '_R2'], '_fastqc.html', True), # path rule 2
        customExpand(['_R1', '_R2'], '.trim.fastq.gz', True), # path rule 4
        customExpand([''], '.bam', True), # path rule 5
        customExpand([''], '.alignlog', True) # path rule 5

################################
# organize sample directories
################################

rule createSampleDirectories:
    input:
        "{sample}_{lane}_{tag}_001.fastq.gz"
    output:
        "{sample}/{sample}_{lane}_{tag}_001.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# QC of fastq files and trim
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
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimPYTHON3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")

################################
# align inserts and fastq screen
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
# sort bam files for filtering
################################

rule sortBam:
    input:
        "{sample}/{sample}_{lane}_001.trim.bam"
    output:
        "{sample}/{sample}_{lane}_001.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={ouput}  SORT_ORDER=coordinate")


################################
# bam filtering with samtools
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
# filter by blacklists
################################

rule filterBamWithBlacklist:
    input:
        "{sample}/{sample}_{lane}_001.trim.st.all.bam"
    output:
        if os.path.exists(config['blacklist']):
            blftBam = "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam"
        filterLog = "{sample}/filtering.log",
        "{sample}/{sample}_{lane}_001.trim.st.all.bam"
    params:
        blacklist = config['blacklist'],
    run:
        # filter by blacklist if provided
        if os.path.exists(config['blacklist']):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > temp.bam") # produces temp file
            shell("samtools view -bh -f 0x2 temp.bam -o {output.blftBam}")
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {output.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {output.filterLog}")

################################
# filter by map quality
################################

# MAKE FUCTIONS FOR CONDITIONAL INPUTS

rule filterBamWithMapQ:
    input:
        inFile = "{sample}/{sample}_{lane}_001.trim.st.all.bam"
        if MAPQ > 0:
            if os.path.exists(config['blacklist']):
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam"
        "{sample}/{sample}_{lane}_001.trim.st.all.bam"
    output:
        if MAPQ > 0:
            if os.path.exists(config['blacklist']):
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.bam",
            else:
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.qft.bam"
        filterLog = "{sample}/filtering.log"
    params:
        blacklist = config['blacklist'],
        mapq = MAPQ
    run:
        if MAPQ > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            shell("samtools view -bh -f 0x2 -q {params.mapq} {input.inFile} -o {output.outFile}")
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {output.filterLog}")
        shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {output.filterLog}")

################################
# histogram with duplicates
################################

# definitely need to break this up more into input/output pairs

rule histogramWithDuplicates_and_removeDuplicates:
    input:
        if MAPQ > 0:
            if os.path.exists(config['blacklist']):
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.bam"
            else:
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.qft.bam"
        else:
            if os.path.exists(config['blacklist']):
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.bam"
            else:
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.bam"
    output:
        if MAPQ > 0:
            if os.path.exists(config['blacklist']):
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.rmdup.bam"
            else:
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.qft.rmdup.bam"
        else:
            if os.path.exists(config['blacklist']):
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.rmdup.bam"
            else:
                outFile = "{sample}/{sample}_{lane}_001.trim.st.all.rmdup.bam"
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log"
    params:
        blacklist = config['blacklist'],
        mapq = 30
    run:
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I={input.inFile} O={output.histDupsLog} H={output.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        shell("picard MarkDuplicates I={input.inFile} O={output.outFile} METRICS_FILE={output.dupsLog} REMOVE_DUPLICATES=true")
        shell("samtools index {output.outFile}")

################################
# bam filtering and enrichment
################################

# definitely need to break this up more into input/output pairs

rule histogramWithNoDuplicates:
    input:
        if MAPQ > 0:
            if os.path.exists(config['blacklist']):
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.qft.rmdup.bam"
            else:
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.qft.rmdup.bam"
        else:
            if os.path.exists(config['blacklist']):
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.blft.rmdup.bam"
            else:
                inFile = "{sample}/{sample}_{lane}_001.trim.st.all.rmdup.bam"
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

