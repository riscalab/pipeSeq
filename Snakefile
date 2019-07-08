#! /bin/bash snakemake
# npagane | 190701 | risca lab | snakefile for ATACseq pipeline

import os
import sys
import datetime

################################
# parameters and functions
################################

# defaults for parameters set in ATACseq.py exectuable file

# determine sample names and sample numbers from the working directory
SAMPLES = []
SAMPLE_NUMS = {}
for base, dirs, files in os.walk("."):
    for fastq in files:
        if fastq.endswith(".fastq.gz") and not fastq.startswith("Undetermined"): # is zipped fastq but NOT undetermined
            tmp = fastq.split(".fastq.gz")[0].split('_')
            SAMPLES.append(tmp[0])
            SAMPLE_NUMS[tmp[0]] = tmp[1]

# generate strucutre of expected files 
def customSeqFileExpand(ext, peak = False): 
    strout = []
    for sample in SAMPLES:
        # check to prepend with sample directory
        if peak:
            primer = sample + '/peakCalls_singles/'
        else: 
            primer = sample + '/'
        # create files
        ftp = primer + sample + '_' + SAMPLE_NUMS[sample] + '_' + config['set'] + ext 
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

################################
# pipeline output check
################################

rule all:
    input:
        customSeqFileExpand(
            conditionalExpand_2(config['mapq'], os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft_summits.bed",
                ".trim.st.all.qft_summits.bed",
                ".trim.st.all.blft_summits.bed",
                ".trim.st.all_summits.bed"
                ),
        True)
    run:
        print('\n###########################')
        print('ATAC-seq pipeline complete')
        print('\n###########################')
        summaryStats()

################################
# run fastq2bam
################################

subworkflow fastq2bam:
    workdir: config["wd"]
    snakefile: "../fastq2bam/Snakefile"

################################
# call peaks to make bed files
################################

rule callPeaks:
    input:
        fastq2bam("{sample}/{sample}_{sampleNum}_{set}.{ext}.rmdup.bam")
    output:
        expand("{{sample}}/peakCalls_singles/{{sample}}_{{sampleNum}}_{{set}}.{{ext}}{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak"])
    params:
        "{sample}/peakCalls_singles/{sample}_{sampleNum}_{set}.{ext}"
    conda:
        "envs/ATACseq.yml" # path relative to snakefile, not working directory
    benchmark:
        "benchmarks/{sample}_{sampleNum}_{set}.{ext}.callPeaks.txt"
    shell:
        "echo 'Calling Peaks...'; macs2 callpeak --nomodel -t {input} -n {params} --nolambda --keep-dup all --call-summits --slocal 10000"

################################
# clean up and summary (7)
################################

# this is executed in the rule all run sequence

def summaryStats():
    with open("ATACseqRunSummary.log", "w") as f: 
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
        f.write('env: fastq2bam\n')
        f.write("SOFTWARE\n")
        f.write("python version: " + str(sys.version_info[0]) + '\n')
        #f.write("pyadapter_trim version: python3 compatible (v1)" + '\n')
        f.write("trim_galore: " + os.popen("trim_galore --v").read() + '\n')
        f.write("fastqc version: " + os.popen("fastqc --version").read() + '\n')
        f.write("bowtie2 version: " + os.popen("bowtie2 --version").read() + '\n')
        f.write("samtools version: " + os.popen("samtools --version").read() + '\n')
        f.write("picard version: 2.20.2-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
        f.write("bedtools version: " + os.popen("bedtools --version").read() + '\n')
        f.write("macs2 version: 2.1.2 <in ATACseq.yml conda env>") # must update whenever ATACseq conda env updated with relevant necessary packages
        f.write("\n\n")
        f.write("PARAMETERS" + '\n')
        f.write("##########\n")
        f.write("genome reference for aligning: " + config["genomeRef"] + '\n')
        f.write("blacklist for filtering: " + config["blacklist"] + '\n')
        f.write("map quality threshold for filtering: " + config["mapq"] + '\n')
        f.write("align command: (bowtie2 -p28 -x {genomeReference} -1 {input.R1} -2 {input.R2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}" + '\n')
        f.write("peak call command: macs2 callpeak --nomodel -t {input} -n {output} --nolambda --keep-dup all --call-summits --slocal 10000")
