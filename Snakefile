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
        #customSeqFileExpand([''],
        #    conditionalExpand_1(os.path.exists(config['blacklist']), 
        #        ".trim.st.all.blft.bam", 
        #        ".trim.st.all.bam"),
        #    True), # rule 5
        #customSeqFileExpand([''],
        #    conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
        #        ".trim.st.all.blft.qft.bam",
        #        ".trim.st.all.qft.bam",
        #        ".trim.st.all.blft.bam",
        #        ".trim.st.all.bam"),
        #    True), # rule 6
        customSeqFileExpand([''],
            conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"),
            True), # rule 7
        customFileExpand('00_source')
    run:
        print('\n###########################')
        print('fastq2bam pipeline complete')
        print('\n###########################')
        shell("mv */.chrM.bam */00_source/")
         mv *.chrM.bam 00_source/
        mv *.all.bam 00_source/
        mv *.st.bam 00_source/
        mv *.rmdup.bam 00_source/
        mv *.st.bam.bai 00_source/


################################
# organize sample directories (0)
################################

rule createSampleDirectories:
    input:
        "{file}.fastq.gz"
    output:
        "{sample}/{file}.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        r1 = "{sample}/{preRun}_R1_001.fastq.gz"
        r2 = "{sample}/{preRun}_R2_001.fastq.gz"
    output:
        "{sample}/{preRun}_R1_001.trim.fastq.gz",
        "{sample}/{preRun}_R2_001.trim.fastq.gz"
    run:
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")
        shell("mv {wildcards.preRun}_R1_001.trim.fastq.gz {wildcards.sample}/")
        shell("mv {wildcards.preRun}_R2_001.trim.fastq.gz {wildcards.sample}/")

################################
# QC of fastq files (2)
################################

rule fastqQC:
    input:
        r1 = "{sample}/{preRun}_R1_001.trim.fastq.gz"
        r2 = "{sample}/{preRun}_R2_001.trim.fastq.gz"
    output:
        "{sample}/{preRun}_R1_001_trim_fastqc.html",
        "{sample}/{preRun}_R2_001_trim_fastqc.html",
        "{sample}/{preRun}_R1_001_trim_fastqc.zip", 
        "{sample}/{preRun}_R2_001_trim_fastqc.zip",
        "{sample}/{preRun}_R1_001.trim.fastq",
        "{sample}/{preRun}_R2_001.trim.fastq"
    run: 
        shell("fastqc -o {wildcards.sample} {input.r1} {input.r2}")
        shell("gunzip {params.r1}")
        shell("gunzip {params.r2}")

################################
# align inserts and fastq screen (2)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = "{sample}/{preRun}_R1_001.trim.fastq",
        unzip2 = "{sample}/{preRun}_R2_001.trim.fastq"
    output:
        "{sample}/{preRun}_R1_001.trim.fastq.gz",
        "{sample}/{preRun}_R2_001.trim.fastq.gz",
        "{sample}/{preRun}_R1_001.trim_screen.html",
        "{sample}/{preRun}_R2_001.trim_screen.html",
        "{sample}/{preRun}_R1_001.trim_screen.png",
        "{sample}/{preRun}_R2_001.trim_screen.png",
        "{sample}/{preRun}_R1_001.trim_screen.txt",
        "{sample}/{preRun}_R2_001.trim_screen.txt",
        bam = "{sample}/{preRun}_001.trim.bam",
        alignLog = "{sample}/{preRun}_001.trim.alignlog"
    params:
        ref = config['genomeRef'],
        screen = "screen.log"
    run:
        shell("(bowtie2 -p28 -x {params.ref} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}")
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {params.screen}")
        shell("mv {params.screen} {wildcards.sample}/")
        shell("mv {wildcards.preRun}_R1_001.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.preRun}_R2_001.trim_screen.* {wildcards.sample}/")
        shell("gzip {input.unzip1}") # zip
        shell("gzip {input.unzip2}") # zip

################################
# sort bam files for filtering (3)
################################

rule sortBam:
    input:
        "{sample}/{sample}_{sampleNu}_001.trim.bam"
    output:
        "{file}.trim.st.bam"
    params:
        "{file}.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")

################################
# bam filtering with samtools (4)
################################

rule filterBam:
    input:
        "{sample}/{file}.trim.st.bam"
    output:
        "{sample}/{file}.trim.st.all.bam"
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
# filter and remove duplicates (5)
################################

rule filter_and_removeDuplicates:
    input:
        "{sample}/{file}.trim.st.all.bam"
    output:
        "{sample}/{file}.{ext}.rmdup.bam"
    params:
        blacklist = config['blacklist'],
        mapq = MAPQ,
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log",
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    run:
        # filter by blacklist if provided
        if os.path.exists(params.blacklist):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > temp.bam") # produces temp file
            ftp = "{wildcards.sample}/{wildcards.file}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 temp.bam -o " + ftp)
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {params.filterLog}")
        else:
            ftp = "{input}"
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")
        # filter by map quality if provided
        if params.mapq > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            tmp = ''.join([ftp.split('.bam')[0] + '.qft' + ftp.split('.bam')[1]])
            shell("samtools view -bh -f 0x2 -q {params.mapq} " + ftp + " -o " + tmp)
            ftp = tmp
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")
        # histogram with duplicates
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        tmp = ''.join([ftp.split('.bam')[0] + '.rmdup' + ftp.split('.bam')[1]])
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={output.histNoDupsLog} H={output.histNoDupsPDF} W=1000 STOP_AFTER=50000")
        # make cleanup directory
        shell("mkdir {wildcards.sample}/00_source")

