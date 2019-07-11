#! /bin/bash snakemake
# npagane | 190701 | risca lab | snakefile rules for fastq2bam pipeline

# include this file to incororporate these rules into a Snakefile for execution

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
def customFileExpand(ext): 
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

# set wildcard constraint (usually 001)
wildcard_constraints:
    set = "\d+"

################################
# align at insertion center (1)
################################

rule ATACoffset:
    input:
        "fastq2bamRunSummary.log", # ensure summary log
        bam = "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.bam"
    output:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bam"
    run:
        shell("bedtools bamtobed -i {input.bam} | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == '+') {$2 = $2 + 4} else if ($6 == '-') {$3 = $3 - 5} print $0}' > {output}")

################################
# call peaks to make bed (2)
################################

rule callPeakSummits:
    input:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bam"
    output:
        expand("{{sample}}/peakCalls_singles/{{sample}}_{{sample_num}}_{{set}}.{{ext}}.atac{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak"]),
        check = "{sample}/{sample}_{sample_num}_{set}.{ext}.atac.CHECK"
    params:
        "{sample}/peakCalls_singles/{sample}_{sample_num}_{set}.{ext}.atac"
    conda:
        "../envs/macs2_python2.yml" # path relative to current file, not working directory
    shell:
        "echo 'Calling Peaks...'; macs2 callpeak --nomodel -t {input} -n {params} --nolambda --keep-dup all --call-summits --slocal 10000; " +
        "touch {output.check}"

################################
# make tracks genomecov (3)
################################

rule bam2bedGraph:
    input:
        bam = "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bam",
        check = "{sample}/{sample}_{sample_num}_{set}.{ext}.atac.CHECK"
    output:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bedgraph"
    run:
        shell("rm {input.check}")
        shell("bedtools genomecov -i {input.bam} -5 -bg > {output}")


################################
# create bigwig for tracks (4)
################################

rule bedGraph2bigWig:
    input:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bedgraph"
    output:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bw"
    params:
        chromSize = config['chromSize']
    run:
        shell("rm {input.check}")
        shell("bedGraphToBigWig {input} {params.chromSize} {output}")

################################
# visualize and analyze tracks (5)
################################

rule visualize_and_analyzeBigWig:
    input:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bw"
    output:
        "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.tab"
    params: 
        cmds = config['igvCommands'],
        bed = "{sample}/{sample}_{sample_num}_{set}.{ext}.rmdup.atac.bedgraph"
    run:
        shell("printf 'load {input}\n" +
             "snapshotDirectory ./{wildcards.sample}/tracks\n';" +
             "cat {params.cmds} > {wildcards.sample}_IGVmds_temp")
        shell("igv -b {wildcards.sample}_IGVcmds_temp")
        shell("rm {wildcards.sample}_cmds_temp") # remove temp file for igv commands
        # statistics
        shell("bigWigAverageOverBed {input} {params.bed} {output}")

################################
# success and summary (6)
################################

rule ATACseqSummary:
    input:
        customFileExpand(
            conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.atac.tab",
                ".trim.st.all.qft.rmdup.atac.tab",
                ".trim.st.all.blft.rmdup.atac.tab",
                ".trim.st.all.rmdup.atac.tab"
            )
        )
    output:
        "ATACseqRunSummary.log"
    run:
        print('\n###########################')
        print('ATAC-seq pipeline complete')
        print('\n###########################')
        with open("{output}", "w") as f: 
            f.write('user: ' + os.environ.get('USER') + '\n')
            f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
            f.write('env: ATACseq\n')
            f.write("SOFTWARE\n")
            f.write("python version: " + str(sys.version_info[0]) + '\n')
            f.write("bedtools version: " + os.popen("bedtools --version").read() + '\n')
            f.write("macs2 version: 2.1.2 <in macs2_python2.yml conda env>\n") # must update if macs2_python2 conda env is updated
            f.write("ucsc tools version: 2 (conda 332)") # must update if new version ever downloaded (shouldnt bc software dependencies)
            f.write("\n\n")
            f.write("PARAMETERS" + '\n')
            f.write("##########\n")
            f.write("chromosome sizes: " + config["chromSize"] + '\n')
            f.write("ivg commands: " + config["ivgCommands"] + '\n')
            f.write("bedGraph command: bedtools genomecov -i {input.bam} -5 -bg > {output}\n")
            f.write("peak call command: macs2 callpeak --nomodel -t {input} -n {output} --nolambda --keep-dup all --call-summits --slocal 10000")
