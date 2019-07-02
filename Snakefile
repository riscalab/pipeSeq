#! /bin/bash snakemake
# npagane | 190701 | risca lab | snakefile for ATACseq pipeline

################################
# parameters and functions
################################

import os
import sys
import datetime

# reassign the file of sample names to a list of the actual samples
SAMPLES = []
with open(config['sampleText']) as ftp:
    lines = ftp.readlines()
    sampleNum = {}
    for line in lines:
        sample = line.strip().split(' ')
        SAMPLES.append(sample[0])
        sampleNum[sample[0]] = str(sample[1])

# defaults for parameters set in ATACseq.py exectuable file

# customized expand function to be compatible with a python dictionary
def customSeqFileExpand(iden, ext, wd = False): 
    strout = []
    for sample in SAMPLES:
        # check to prepend with sample directory
        if wd:
            primer = sample + '/peakCalls_singles/'
        else:
            primer = ''
        # create files
        for i in iden:
            ftp = primer + sample + '_S' + sampleNum[sample] + i + '_' + config['set'] + ext 
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
        customSeqFileExpand([''],
            conditionalExpand_2(config['mapq'], os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft_summits.bed",
                ".trim.st.all.qft_summits.bed",
                ".trim.st.all.blft_summits.bed",
                ".trim.st.all_summits.bed"),
            True), # rule 7
    run:
        print('\n###########################')
        print('ATAC-seq pipeline complete')
        print('\n###########################')
        summaryStats()

################################
# run fastq2bam
################################

subworkflow otherworkflow:
    snakefile: "../fastq2bam/Snakefile"

################################
# call peaks to make bed files
################################

rule callPeaks:
    input:
        otherworkflow("{sample}/{sample}_{sampleNum}_{set}.{ext}.rmdup.bam")
    output:
        expand("{{sample}}/peakCalls_singles/{{sample}}_{{sampleNum}}_{{set}}.{{ext}}{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak"])
    params:
        "{sample}/peakCalls_singles/{sample}_{sampleNum}_{set}.{ext}"
    conda:
        "envs/ATACseq.yml"
    shell:
        "echo 'Calling Peaks...'; macs2 callpeak --nomodel -t {input} -n {params} --nolambda --keep-dup all --call-summits --slocal 10000"

################################
# clean up and summary (7)
################################

# this is executed in the rule all run sequence

def summaryStats():
    with open("ATACseqSummary.txt", "w") as f: 
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
        f.write("SOFTWARE\n")
        f.write("########\n")
        f.write("bedtools version: " + os.popen("bedtools --version").read() + '\n')
        f.write("\n\n")
        f.write("PARAMETERS" + '\n')
        f.write("##########\n")
        f.write("genome reference for aligning: " + config["genomeRef"] + '\n')
        f.write("blacklist for filtering: " + config["blacklist"] + '\n')
        f.write("map quality threshold for filtering: " + config["mapq"] + '\n')
