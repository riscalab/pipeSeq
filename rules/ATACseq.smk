#! usr/bin/env snakemake
## npagane | 190701 | risca lab | ATACseq pipeline rules

# include this file to incororporate these rules into a Snakefile for execution

import os
import sys
import datetime

################################
# parameters and functions
################################

# defaults for parameters set in ATACseq.py exectuable file

# functions are found in helper.py file

# determine if an igvCommands file was passed
if not os.path.exists(config['igvCommands']):
    config['igvCommands'] = workflow.basedir + "/scripts/igvCommands"

wildcard_constraints:
    post_tag = "\d+"

################################
# align at insertion center (1)
################################

rule ATACoffset:
    input:
        "fastq2bamRunSummary.log", # ensure summary log
        bam = "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bam"
    params:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.temp.bam"
    threads: 2
    run:
        shell("samtools index {input.bam}") # suppress the pysam/htslib warning about the index file
        shell("alignmentSieve --numberOfProcessors {threads} --ATACshift --bam {input.bam} -o {params}")
        shell("samtools sort -O bam -o {output} {params}") #sort (dont use picard it is too strict about sam formatting)
        shell("samtools index {output}") # regenerate index file
        shell("rm {params}")

################################
# call peaks to make bed (2)
################################

rule callPeakSummits:
    input:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bam"
    output:
        expand("{{sample}}/peakCalls_singles/{{pre_tag}}_{{post_tag}}.{{ext}}.atac{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak"]),
        check = "{sample}/{pre_tag}_{post_tag}.{ext}.atac.CHECK"
    params:
        "{sample}/peakCalls_singles/{pre_tag}_{post_tag}.{ext}.atac"
    conda:
        "../envs/macs2_python2.yml" # path relative to current file, not working directory
    shell:
        "echo 'Calling Peaks...'; macs2 callpeak --nomodel -f BAM -t {input} -n {params} --nolambda --keep-dup all --call-summits --slocal 10000; " +
        "touch {output.check}"

################################
# make tracks genomecov (3)
################################

rule bam2bedGraph:
    input:
        bam = "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bam",
        check = "{sample}/{pre_tag}_{post_tag}.{ext}.atac.CHECK"
    output:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bedgraph"
    params:
        chromSize = config['chromSize']
    run:
        shell("rm {input.check}")
        shell("bedtools genomecov -ibam {input.bam} -5 -bg > {output}")


################################
# create bigwig for tracks (4)
################################

rule bedGraph2bigWig:
    input:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bedgraph"
    output:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bw"
    params:
        chromSize = config['chromSize'],
        st = "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.st.bedgraph"
    run:
        shell("LC_COLLATE=C sort -k1,1 -k2,2n {input} > {params.st}")
        shell("bedGraphToBigWig {params.st} {params.chromSize} {output}")

################################
# visualize and analyze tracks (5)
################################

rule visualize_and_analyzeBigWig:
    input:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bw",
    output:
        "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.tab"
    params: 
        cmds = config['igvCommands'],
        bed = "{sample}/{pre_tag}_{post_tag}.{ext}.rmdup.atac.bedgraph",
        tmp = "{sample}_IGVcmds_tmp"
    run:
        shell("printf 'load {input}\n" + # load genome around here
             "snapshotDirectory {wildcards.sample}/tracks\n';" +
             "cat {params.cmds} > {params.tmp}")
        shell("igv -b {params.tmp}")
        shell("rm {params.tmp}") # remove temp file for igv commands
        # statistics
        shell("bigWigAverageOverBed {input} {params.bed} {output}")

################################
# success and summary (6)
################################

rule ATACseqSummary:
    input:
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.atac.tab",
                ".trim.st.all.qft.rmdup.atac.tab",
                ".trim.st.all.blft.rmdup.atac.tab",
                ".trim.st.all.rmdup.atac.tab"
            ), config['exclude']
        )
    output:
        "ATACseqRunSummary.log"
    run:
        print('\n###########################')
        print('ATAC-seq pipeline complete')
        print('\n###########################')
        with open(output[0], "w") as f: 
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
            f.write("ivg commands: " + config["igvCommands"] + '\n')
            f.write("bedGraph command: bedtools genomecov -i {input.bam} -5 -bg > {output}\n")
            f.write("peak call command: macs2 callpeak --nomodel -t {input} -n {output} --nolambda --keep-dup all --call-summits --slocal 10000")
            f.write("SUMMARY\n")
            f.write("#######\n")
