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
    config['blacklist'] = "/rugpfs/fs0/risc_lab/scratch/nvelez/blacklists/ATAC_blacklist.bed"
    #config['blacklist'] = ""

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

# set wildcard constraint (usually 001)
wildcard_constraints:
    set = "\d+"

# map quality
MAPQ = 30 # maybe make this a prompt at somepoint.. but for now, this is fine

################################
# pipeline output check
################################

rule all:
    input:
        customSeqFileExpand(TAGS, '.fastq.gz'),
        customSeqFileExpand([''],
            conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"),
            True), # rule 7
    run:
        print('\n###########################')
        print('fastq2bam pipeline complete')
        print('\n###########################')
        # will sadly have to loop through to clean instead of using fun asterisks
        for s in SAMPLES:
            #shell("mv " + s +"/*.chrM.bam " + s + "/00_source/") # this usually does not exist
            shell("mv " + s +"/*.all.bam " + s + "/00_source/")
            shell("mv " + s +"/*.st.bam.bai " + s + "/00_source/") # am not moving rmdup files

################################
# organize sample directories (0)
################################

rule createSampleDirectories:
    input:
        "{sample}_{sampleNum}_{tag}_{set}.fastq.gz"
    output:
        "{sample}/{sample}_{sampleNum}_{tag}_{set}.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        r1 = "{sample}/{sample}_{sampleNum}_R1_{set}.fastq.gz",
        r2 = "{sample}/{sample}_{sampleNum}_R2_{set}.fastq.gz"
    output:
        expand("{{sample}}/{{sample}}_{{sampleNum}}_{run}_{{set}}.trim.fastq.gz", run=["R1", "R2"])
    params:
        r1 = "{sample}_{sampleNum}_R1_{set}.trim.fastq.gz",
        r2 = "{sample}_{sampleNum}_R2_{set}.trim.fastq.gz"
    run:
        shell("/rugpfs/fs0/risc_lab/store/risc_soft/pyadapter_trim/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")
        shell("mv {params.r1} {wildcards.sample}/")
        shell("mv {params.r2} {wildcards.sample}/")

################################
# QC of fastq files (2)
################################

rule fastqQC:
    input:
        r1 = "{sample}/{sample}_{sampleNum}_R1_{set}.trim.fastq.gz",
        r2 = "{sample}/{sample}_{sampleNum}_R2_{set}.trim.fastq.gz"
    output:
        expand("{{sample}}/{{sample}}_{{sampleNum}}_{run}_{{set}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", ".fastqc.zip", ".fastq"])
    params:
        r1 = "{sample}/{sample}_{sampleNum}_R1_{set}.trim.fastq.gz",
        r2 = "{sample}/{sample}_{sampleNum}_R2_{set}.trim.fastq.gz"
    run: 
        shell("fastqc -o {wildcards.sample} {input.r1} {input.r2}")
        shell("gunzip {params.r1}")
        shell("gunzip {params.r2}")

################################
# align inserts and fastq screen (2)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = "{sample}/{sample}_{sampleNum}_R1_{set}.trim.fastq",
        unzip2 = "{sample}/{sample}_{sampleNum}_R2_{set}.trim.fastq"
    output:
        expand("{{sample}}/{{sample}}_{{sampleNum}}_{run}_{{set}}.trim{end}", run=["R1", "R2"], end=["_screen.html", "_screen.png", "_screen.txt"]),
        bam = "{sample}/{sample}_{sampleNum}_{set}.trim.bam",
        alignLog = "{sample}/{sample}_{sampleNum}_{set}.trim.alignlog"
    params:
        ref = config['genomeRef'],
        screen = "screen.log", 
    run:
        shell("(bowtie2 -p28 -x {params.ref} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}")
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {params.screen}")
        shell("mv {params.screen} {wildcards.sample}/")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R1_{wildcards.set}.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R2_{wildcards.set}.trim_screen.* {wildcards.sample}/")
        shell("gzip {input.unzip1}") # zip
        shell("gzip {input.unzip2}") # zip

################################
# sort bam files for filtering (3)
################################

rule sortBam:
    input:
        "{sample}/{sample}_{sampleNum}_{set}.trim.bam"
    output:
        "{sample}/{sample}_{sampleNum}_{set}.trim.st.bam"
    params:
        "{sample}_{sampleNum}_{set}.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")

################################
# bam filtering with samtools (4)
################################

rule filterBam:
    input:
        "{sample}/{sample}_{sampleNum}_{set}.trim.st.bam"
    output:
        "{sample}/{sample}_{sampleNum}_{set}.trim.st.all.bam"
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
        "{sample}/{sample}_{sampleNum}_{set}.trim.st.all.bam"
    output:
        "{sample}/{sample}_{sampleNum}_{set}.{ext}.rmdup.bam"
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
            ftp = "{wildcards.sample}/{wildcards.sample}_{wildcards.sampleNum}_{wildcards.set}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 temp.bam -o " + ftp)
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
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        tmp = ftp.split('.bam')[0] + '.rmdup.bam'
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histNoDupsLog} H={params.histNoDupsPDF} W=1000 STOP_AFTER=50000")
        # make cleanup directory
        shell("mkdir {wildcards.sample}/00_source")

