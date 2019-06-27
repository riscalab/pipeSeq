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
    cnt = 1
    for line in lines:
        sample = line.strip()
        SAMPLES.append(sample)
        sampleNum[sample] = str(cnt)
        cnt += 1

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
        customSeqFileExpand(TAGS, '.fastq.gz', True), # rule 1
        customSeqFileExpand(['_R1', '_R2'], '.trim.fastq.gz', True), # rule 2
        customSeqFileExpand([''], '.trim.bam', True), # rule 3
        customSeqFileExpand([''], '.trim.st.bam', True), # rule 4
        customSeqFileExpand([''], '.trim.st.all.bam', True), # rule 5
        customSeqFileExpand([''],
            conditionalExpand_1(os.path.exists(config['blacklist']), 
                ".trim.st.all.blft.bam", 
                ".trim.st.all.bam"),
            True), # rule 6
        customSeqFileExpand([''],
            conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.bam",
                ".trim.st.all.qft.bam",
                ".trim.st.all.blft.bam",
                ".trim.st.all.bam"),
            True), # rule 7
        #customSeqFileExpand([''],
        #    conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
        #        ".trim.st.all.blft.qft.rmdup.bam",
        #        ".trim.st.all.qft.rmdup.bam",
        #        ".trim.st.all.blft.rmdup.bam",
        #        ".trim.st.all.rmdup.bam"),
        #    True), # rule 8
        #customFileExpand("hist_datawithoutdups.log") # rule 9

################################
# organize sample directories (1)
################################

rule createSampleDirectories:
    input:
        "{sample}_{sampleNum}_{tag}_001.fastq.gz"
    output:
        "{sample}/{sample}_{sampleNum}_{tag}_001.fastq.gz"
    shell:
        "mv {input} {output}"


################################
# QC of fastq files and trim (2)
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
# align inserts and fastq screen (3)
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
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {params.screen} | mv {params.screen} {wildcards.sample}/") # OK LOOK HERE fastqc screen is weird... should i just add all genomes 
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R1_001.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.sample}_{wildcards.sampleNum}_R2_001.trim_screen.* {wildcards.sample}/")
        shell("gzip {input.unzip1}") # zip
        shell("gzip {input.unzip2}") # zip

################################
# sort bam files for filtering (4)
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
# bam filtering with samtools (5)
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
        shell("echo 'Filtering file {input} by script sh02a_filter_bam.sh' >> {params.filterLog} | mv {params.filterLog} {wildcards.sample}/")

################################
# filter by blacklists (6)
################################

rule filterBamWithBlacklist:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all.bam"
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}.bam",
        blft = "{sample}/{sample}_{sampleNum}_blftStep{blft,.*}.txt" # hardcode in input/output sequence in case of no filter
    params:
        blacklist = config['blacklist'],
        filterLog = "{sample}/filtering.log"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", "blft"]])
    run:
        # filter by blacklist if provided
        if os.path.exists(params.blacklist):
            shell("echo 'sh02a_filter_bam.sh: Removing blacklisted reads")
            shell("bedtools intersect -v -abam {input} -b {params.blacklist} -wa > temp.bam") # produces temp file
            shell("samtools view -bh -f 0x2 temp.bam -o {wildcards.sample}/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all.blft.bam")
            shell("echo 'Blacklist filtered using file {params.blacklist}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")
        shell("touch {output.blft}")

################################
# filter by map quality (7)
################################

rule filterBamWithMapQ:
    input:
        inFile = "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}.bam",
        blft = "{sample}/{sample}_{sampleNum}_blftStep{blft,.*}.txt" # hardcode in input/output sequence in case of no filter
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}{qft,.*}.bam",
        mapq = "{sample}/{sample}_{sampleNum}_mapqStep{blft,.*}{qft,.*}.txt" # hardcode in input/output sequence in case of no filter
    params:
        mapq = MAPQ,
        filterLog = "{sample}/filtering.log"
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", "blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", "qft"]])
    run:
        print('hi')
        if params.mapq > 0:
            shell("echo 'sh02a_filter_bam.sh: Removing low quality reads'")
            shell("samtools view -bh -f 0x2 -q {params.mapq} {input.inFile} -o {wildcards.sample}/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all{wildcards.blft}.qft.bam")
            shell("echo 'Filtered with mapping quality filter {params.mapq}.' >> {params.filterLog}")
        else:
            shell("echo 'sh02a_filter_bam.sh: Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")
        shell("rm {input.blft}")
        shell("touch {output.mapq}")

################################
# histogram with duplicates (8)
################################
'''
rule histogramWithDuplicates_and_removeDuplicates:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}{qft,.*}.bam"
        #inFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.qft.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.qft.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.bam")
    output:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}{qft,.*}.rmdup.bam",
        #outFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.qft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.qft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.rmdup.bam"),
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log"
    params:
        blacklist = config['blacklist'],
        mapq = MAPQ
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", "blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", "qft"]])
    run:
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetricCs I={input} O={output.histDupsLog} H={output.histDupsPDF} W=1000 STOP_AFTER=50000")
        # remove duplicates
        shell("echo 'sh02b_remove_dups_estimate_diversity.sh: Removing duplicates'")
        shell("picard MarkDuplicates I={input} O={wildcards.sample}/{wildcards.sample}_{wildcards.sampleNum}_001.trim.st.all{wildcards.blft}{wildcards.qft}.rmdup.bam METRICS_FILE={output.dupsLog} REMOVE_DUPLICATES=true")
        #shell("samtools index {output.outFile}") # moved this down a rule bc it should operate on input not output

################################
# histogram with removed dups (9)
################################

rule histogramWithNoDuplicates:
    input:
        "{sample}/{sample}_{sampleNum}_001.trim.st.all{blft,.*}{qft,.*}.rmdup.bam"
        #inFile = conditionalExpand_2(MAPQ, os.path.exists(config['blacklist']),
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.qft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.qft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.blft.rmdup.bam",
        #    "{sample}/{sample}_{sampleNum}_001.trim.st.all.rmdup.bam")
    output:
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    params:
        blacklist = config['blacklist'],
        mapq = 30
    wildcard_constraints:
        blft = '|'.join([re.escape(x) for x in ["", "blft"]]),
        qft = '|'.join([re.escape(x) for x in ["", "qft"]])
    run:
        shell("samtools index {input}")
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I={input} O={output.histNoDupsLog} H={output.histNoDupsPDF} W=1000 STOP_AFTER=50000")


################################
# clean up ??
################################

rule cleanUp:
    input:
         "{sample}/hist_data_withoutdups.log"
    output:
        "{sample}/00_source/"
'''
