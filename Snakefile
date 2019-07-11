#! /bin/bash snakemake
# npagane | 190618 | risca lab | snakefile for fastq2bam pipeline

import os
import sys
import datetime

################################
# parameters and functions
################################

# defaults for parameters set in fastq2bam.py exectuable file

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
def customSeqFileExpand(ext): 
    strout = []
    for sample in SAMPLES:
        ftp = sample + '/' + sample + '_' + SAMPLE_NUMS[sample] + '_' + config['set'] + ext 
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

################################
# pipeline output check
################################

rule all:
    input:
        customSeqFileExpand(
            conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"
            )
        )
    run:
        print('\n###########################')
        print('fastq2bam pipeline complete')
        print('\n###########################')
        summaryStats()

################################
# rules for fastq2bam
################################

include: "./rules/fastq2bam.smk"

################################
# clean up and summary 
################################

def summaryStats():
    with open("fastq2bamRunSummary.log", "w") as f:
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
        f.write("PARAMETERS" + '\n')
        f.write("##########\n")
        f.write("genome reference for aligning: " + config["genomeRef"] + '\n')
        f.write("blacklist for filtering: " + config["blacklist"] + '\n')
        f.write("map quality threshold for filtering: " + config["mapq"] + '\n')
        f.write("align command: (bowtie2 -p28 -x {genomeReference} -1 {input.R1} -2 {input.R2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}" + '\n')
