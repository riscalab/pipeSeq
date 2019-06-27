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
    sampleNum = {}
    for line in lines:
        sample = line.strip().split(' ')
        SAMPLES.append(sample[0])
        sampleNum[sample[0]] = str(sample[1])

# defaults for parameters
if not os.path.exists(config['genomeRef']):
    config['genomeRef'] = "/rugpfs/fs0/risc_lab/scratch/nvelez/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"

if not os.path.exists(config['blacklist']):
    #config['blacklist'] = "/rugpfs/fs0/risc_lab/scratch/nvelez/blacklists/ATAC_blacklist.bed"
    config['blacklist'] = ""

if not os.path.exists(config['fastaRef']):
    config['fastaRef'] = "/rugpfs/fs0/risc_lab/scratch/nvelez/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa"


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
            ftp = primer + sample + '_S' + sampleNum[sample] + i + '_001' + ext 
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
        customSeqFileExpand(TAGS, '.fastq.gz', True), # rule 0
        customSeqFileExpand(['_R1', '_R2'], '.trim.fastq.gz', True), # rule 1
        customSeqFileExpand([''], '.trim.bam', True), # rule 2
        customSeqFileExpand([''], '.trim.st.bam', True), # rule 3
        customSeqFileExpand([''], '.trim.st.all.bam', True), # rule 4
        customSeqFileExpand([''],
            conditionalExpand_1(os.path.exists(config['blacklist']), 
                ".trim.st.all.blft.bam", 
                ".trim.st.all.bam"),
            True), # rule 5
        customSeqFileExpand([''],
            conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.bam",
                ".trim.st.all.qft.bam",
                ".trim.st.all.blft.bam",
                ".trim.st.all.bam"),
            True), # rule 6
        customSeqFileExpand([''],
            conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"),
            True), # rule 7
        customFileExpand('00_source')

################################
# organize sample directories (0)
################################

rule createSampleDirectories:
    input:
        "{sample}_{sampleNum}_{tag}_001.fastq.gz"
    output:
        "{sample}/{sample}_{sampleNum}_{tag}_001.fastq.gz"
    shell:
        "mv {input} {output}"


################################
# QC of fastq files and trim (1)
################################

rule fastqQC_and_trimAdapters:
    input:
        r1 = "{sample}/{sample}_{sampleNum}_R1_001.fastq.gz",
        r2 = "{sample}/{sample}_{sampleNum}_R2_001.fastq.gz"
    output:
        "{sample}/{sample}_{sampleNum}_R1_001_fastqc.html",
        "{sample}/{sample}_{sampleNum}_R2_001_fastqc.html",
        "{sample}/{sample}_{sampleNum}_R1_001_fastqc.zip", 
        "{sample}/{sample}_{sampleNum}_R2_001_fastqc.zip",
        "{sample}/{sample}_{sampleNum}_R1_001.trim.fastq",
        "{sample}/{sample}_{sampleNum}_R2_001.trim.fastq"
    params:
        r1 = "{sample}/{sample}_{sampleNum}_R1_001.trim.fastq.gz",
        r2 = "{sample}/{sample}_{sampleNum}_R2_001.trim.fastq.gz"
    run: 
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")
        shell("fastqc -o {wildcards.sample} {input.r1} {input.r2}")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R1_001.trim.fastq.gz {params.r1}")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R2_001.trim.fastq.gz {params.r2}")
        shell("gunzip {params.r1}")
        shell("gunzip {params.r2}")

################################
# align inserts and fastq screen (2)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = "{sample}/{sample}_{sampleNum}_R1_001.trim.fastq",
        unzip2 = "{sample}/{sample}_{sampleNum}_R2_001.trim.fastq"
    output:
        "{sample}/{sample}_{sampleNum}_R1_001.trim.fastq.gz",
        "{sample}/{sample}_{sampleNum}_R2_001.trim.fastq.gz",
        "{sample}/{sample}_{sampleNum}_R1_001.trim_screen.html",
        "{sample}/{sample}_{sampleNum}_R2_001.trim_screen.html",
        "{sample}/{sample}_{sampleNum}_R1_001.trim_screen.png",
        "{sample}/{sample}_{sampleNum}_R2_001.trim_screen.png",
        "{sample}/{sample}_{sampleNum}_R1_001.trim_screen.txt",
        "{sample}/{sample}_{sampleNum}_R2_001.trim_screen.txt",
        bam = "{sample}/{sample}_{sampleNum}_001.trim.bam",
        alignLog = "{sample}/{sample}_{sampleNum}_001.trim.alignlog"
    params:
        ref = config['genomeRef'],
        screen = "screen.log"
    run:
        shell("(bowtie2 -p28 -x {params.ref} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}")
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {params.screen}")
        shell("mv {params.screen} {wildcards.sample}/")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R1_001.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R2_001.trim_screen.* {wildcards.sample}/")
        shell("gzip {input.unzip1}") # zip
        shell("gzip {input.unzip2}") # zip

################################
# sort bam files for filtering (3)
################################

rule sortBam:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.bam"
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.bam"
    params:
        "{sample}_{sampleNum}_001.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")

################################
# bam filtering with samtools (4)
################################

rule filterBam:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.bam"
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all.bam"
    params:
        chrs = "samtools view -H *.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v _gl | grep -v Y | grep -v hap | grep -v random | grep -v v1 | grep -v v2",
        filterLog = "filtering.log"
    run:
        shell("samtools index {input}")
        shell("echo 'sh02a_filter_bam.sh: Removing reads from unwanted chromosomes and scaffolds'")
        shell("samtools view -b {input} `echo {params.chrs}` > {output}")
        shell("echo 'Filtering file {input} by script sh02a_filter_bam.sh' >> {params.filterLog}")
        shell("mv {params.filterLog} {wildcards.sample}/")

################################
# filter by blacklists (5)
################################

rule filterBamWithBlacklist:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all.bam"
    output:
        "{sample}/BF/{sample}_{sampleNum}_001.trim.st.all{blft}.bam"
    params:
        blacklist = config['blacklist'],
        filterLog = "{sample}/filtering.log"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", ".blft"]])
    run:
        # filter by blacklist if provided
        if os.path.exists(params.blacklist):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > temp.bam") # produces temp file
            shell("samtools view -bh -f 0x2 temp.bam -o {wildcards.sample}/BF/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all.blft.bam")
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")

################################
# filter by map quality (6)
################################

rule filterBamWithMapQ:
    input:
        "{sample}/BF/{sample}_{sampleNum}_001.trim.st.all{blft}.bam"
    output:
        "{sample}/BFQF{sample}_{sampleNum}_001.trim.st.all{blft}{qft}.bam"
    params:
        mapq = MAPQ,
        filterLog = "{sample}/filtering.log"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", ".blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", ".qft"]])
    run:
        if params.mapq > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            shell("samtools view -bh -f 0x2 -q {params.mapq} {input} -o {wildcards.sample}/BFQF/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all{wildcards.blft}.qft.bam")
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")

################################
# histogram with duplicates (7)
################################

rule histogramWithDuplicates_and_removeDuplicates:
    input:
        "{sample}/BFQF{sample}_{sampleNum}_001.trim.st.all{blft}{qft}.bam"
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft}{qft}.rmdup.bam"
    params:
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", ".blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", ".qft"]])
    run:
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetricCs I={input} O={params.histDupsLog} H={params.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        shell("picard MarkDuplicates I={input} O={wildcards.sample}/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all{wildcards.blft}{wildcards.qft}.rmdup.bam METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        #shell("samtools index {output.outFile}") # moved this down a rule bc it should operate on input not output

################################
# histogram with removed dups (8)
################################

rule histogramWithNoDuplicates:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft}{qft}.rmdup.bam"
    output:
        "{sample}/{samle}_{sampleNum}{blft}{qft}DONE.txt" # MAKE COMPREHENSIVE PLOTS HERE
    params:
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", ".blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", ".qft"]])
    run:
        shell("samtools index {input}")
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I={input} O={output.histNoDupsLog} H={output.histNoDupsPDF} W=1000 STOP_AFTER=50000")


################################
# clean up ??
################################

rule cleanUp:
    input:
        "{sample}/*DONE.txt",
        bf = "{sample}/BF/{file}.bam",
        bfqf = "{sample}/BFQF/{file}.bam"
    output:
        "{sample}/{file}.bam"
    wildcard_constraints:
        sample = "[A-Za-z0-9]+"
    run:
        shell("mv {input.bf} {ouput}")
        shell("mv {input.bfqf} {output}")
        shell("mkdir {wildcards.sample} 00_source")

################################
# success and summary
################################
'''
rule successSummary:
    input:
    output:
        "{sample}/{file}.bam"
    run:
        shell("mv {input.bf} {ouput}")
        shell("mv {input.bfqf} {output}")
        shell("mkdir {wildcards.sample} 00_source")
'''
