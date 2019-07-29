#! usr/bin/env snakemake
## npagane | 190701 | risca lab | fastq2bam pipeline rules

# include this file to incororporate these rules into a Snakefile for execution

import os
import numpy as np
import sys
import helper
import datetime

################################
# parameters and functions
################################

# defaults for parameters set in fastq2bam.py exectuable file

# functions are found in helper.py file

# determine if there are index fastq files
if config['index'] == 'True':
    TAGS = ['R1', 'R2', 'I1', 'I2'] 
else:
    TAGS = ['R1', 'R2']

wildcard_constraints:
    post_tag = "\d+"

################################
# commands with custom flags
################################      

align = "(bowtie2 -p{threads} -x {config[genomeRef]} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS -o {output.bam}) 2>{output.alignLog}"

################################
# organize sample directories (0)
################################

rule createSampleDirectories:
    input:
        "{pre_tag}_{tag}_{post_tag}.fastq.gz"
    output:
        "{sample}/{pre_tag}_{tag}_{post_tag}.fastq.gz"
    shell:
        "mv {input} {output}"

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        expand("{{sample}}/{{pre_tag}}_{tag}_{{post_tag}}.fastq.gz", tag=TAGS), # move all files to sample dirs
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.fastq.gz"
    output:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    params:
        r1 = "{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    run:
        shell(workflow.basedir + "/scripts/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {wildcards.sample}/adapter_trim.log")
        shell("mv {params.r1} {wildcards.sample}/")
        shell("mv {params.r2} {wildcards.sample}/")

################################
# QC of fastq files (2)
################################

rule fastqc:
    input:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    output:
        expand("{{sample}}/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", "_fastqc.zip", ".fastq"])
    params:
        r1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    threads: 2
    run:
        shell("fastqc -t {threads} -o {wildcards.sample} {input.r1} {input.r2}")
        shell("unpigz {params.r1} {params.r2}")

################################
# align inserts and fastq screen (3)
################################

rule alignInserts_and_fastqScreen:
    input:
        unzip1 = "{sample}/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = "{sample}/{pre_tag}_R2_{post_tag}.trim.fastq",
    output:
        expand("{{sample}}/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_screen.html", "_screen.png", "_screen.txt"]),
        bam = "{sample}/{pre_tag}_{post_tag}.trim.bam",
        alignLog = "{sample}/{pre_tag}_{post_tag}.trim.alignlog"
    params:
        screen = "screen.log"
    threads: 8
    run:
        shell(align) # align command defined above
        shell("fastq_screen --aligner bowtie2 {input.unzip1} {input.unzip2}  > {wildcards.sample}_{params.screen}")
        shell("mv {wildcards.sample}_{params.screen} {wildcards.sample}/{params.screen}")
        shell("mv {wildcards.pre_tag}_R1_{wildcards.post_tag}.trim_screen.* {wildcards.sample}/")
        shell("mv {wildcards.pre_tag}_R2_{wildcards.post_tag}.trim_screen.* {wildcards.sample}/")
        shell("pigz {input.unzip1} {input.unzip2}") # zip

################################
# sort bam files for filtering (4)
################################

rule sortBam:
    input:
        "{sample}/{pre_tag}_{post_tag}.trim.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.trim.st.bam"
    params:
        "{pre_tag}_{post_tag}.trim.st.bam"
    run:
        shell("picard SortSam  I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")

################################
# bam filtering with samtools (5)
################################

rule filterBam:
    input:
        "{sample}/{pre_tag}_{post_tag}.trim.st.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}.trim.st.all.bam"
    params:
        chrs = "samtools view -H {sample}/{pre_tag}_{post_tag}.trim.st.bam | grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM | grep -v _gl | grep -v Y | grep -v hap | grep -v random | grep -v v1 | grep -v v2",
        chrM = "{sample}/{pre_tag}_{post_tag}.trim.st.chrM.bam",
        filterLog = "filtering.log"
    run:
        shell("samtools index {input}")
        shell("echo 'Removing reads from unwanted chromosomes and scaffolds'")
        shell("samtools view -b {input} `echo {params.chrs}` > {output}")
        shell("samtools view -b {input} chrM > {params.chrM}")
        shell("echo 'Filtering file {input} by rules filterBam and filter_removeDups_and_enrichTSS' >> {wildcards.sample}_{params.filterLog}")
        shell("mv {wildcards.sample}_{params.filterLog} {wildcards.sample}/{params.filterLog}")

################################
# blacklist and map quality filter, 
# remove duplicates, and TSS enrichment (6)
################################

rule filter_removeDups_and_enrichTSS:
    input:
        "{sample}/{pre_tag}_{post_tag}.trim.st.all.bam"
    output:
        "{sample}/{pre_tag}_{post_tag}{ext,.*}.rmdup.bam"
    params:
        filterLog = "{sample}/filtering.log",
        histDupsLog = "{sample}/hist_data_withdups.log",
        histDupsPDF = "{sample}/hist_graphwithdups.pdf",
        dupsLog = "{sample}/dups.log",
        histNoDupsLog = "{sample}/hist_data_withoutdups.log",
        histNoDupsPDF = "{sample}/hist_graphwithoutdups.pdf"
    run:
        # filter by blacklist if provided
        if os.path.exists(config["blacklist"]):
            shell("echo 'Removing blacklisted reads'")
            shell("bedtools intersect -v -abam {input} -b {config[blacklist]} -wa > {wildcards.sample}_temp.bam") # produces temp file
            ftp = "{wildcards.sample}/{wildcards.pre_tag}_{wildcards.post_tag}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 {wildcards.sample}_temp.bam -o " + ftp)
            shell("rm {wildcards.sample}_temp.bam") # remove temp file
            shell("echo 'Blacklist filtered using file {config[blacklist]}.' >> {params.filterLog}")
        else:
            ftp = input
            shell("echo 'Blacklist file not found or specified. Skipping blacklist filter.'")
            shell("echo 'Did not filter by blacklist.' >> {params.filterLog}")
        # filter by map quality if provided
        if int(config["mapq"]) > 0:
            shell("echo 'Removing low quality reads'")
            tmp = str(ftp).split('.bam')[0] + '.qft.bam'
            shell("samtools view -bh -f 0x2 -q {config[mapq]} " + str(ftp) + " -o " + tmp)
            ftp = tmp
            shell("echo 'Filtered with mapping quality filter {config[mapq]}.' >> {params.filterLog}")
        else:
            shell("echo 'Mapping quality threshold not specified, quality filter skipped'")
            shell("echo 'Did not filter by mapping quality.' >> {params.filterLog}")
        # histogram with duplicates
        shell("echo 'Histogram with duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=600 STOP_AFTER=5000000")
        # remove duplicates
        shell("echo 'Removing duplicates'")
        tmp = ftp.split('.bam')[0] + '.rmdup.bam'
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histNoDupsLog} H={params.histNoDupsPDF} W=600 STOP_AFTER=5000000")
        # TSS enrichment if provided
        if os.path.exists(config["TSS"]):
            tmp = ftp + '.RefSeqTSS'
            shell(workflow.basedir + "/scripts/pyMakeVplot_css_v01.py -a " + ftp + " -b {config[TSS]} -e 2000 -p ends -s 5 -v -u --atac -o " + tmp)
        else:
            shell("echo 'TSS BED file not provided. not creating TSS enrichment profile'")
        # cleanup directory
        shell("if [ ! -d {wildcards.sample}/00_source ]; then mkdir {wildcards.sample}/00_source; fi") 
        shell("mv {wildcards.sample}/*.all.bam {wildcards.sample}/00_source/")
        shell("mv {wildcards.sample}/*.chrM.bam {wildcards.sample}/00_source/")

################################
# combine insert plots and summary (7)
################################

rule fastq2bamSummary:
    input:
        helper.customFileExpand(
            helper.conditionalExpand_2(int(config['mapq']), os.path.exists(config['blacklist']),
                ".trim.st.all.blft.qft.rmdup.bam",
                ".trim.st.all.qft.rmdup.bam",
                ".trim.st.all.blft.rmdup.bam",
                ".trim.st.all.rmdup.bam"
            ), config['exclude']
        )
    output:
        "fastq2bamRunSummary.log"
    params:
        files = np.unique(list(helper.findFiles(config['exclude']).keys())), # order the samples
        temp = "tempSummary_fastq.log"
    run:
        # make nice pdf of insert distributions
        libs = "libsTEMP_fastq.txt"
        with open(libs, "w") as f:
            for i in params.files:
                f.write(i + '\n')
        shell("Rscript " + workflow.basedir + "/scripts/plotisds_v2.R " + libs + " hist_data_withoutdups")
        shell("rm " + libs)
        print('\n###########################')
        print('fastq2bam pipeline complete')
        print('\n###########################')
        with open(output[0], "w") as f:
            f.write('user: ' + os.environ.get('USER') + '\n')
            f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
            f.write("SOFTWARE\n")
            f.write("########\n")
            f.write("python version: " + str(sys.version_info[0]) + '\n')
            f.write("pyadapter_trim version: python3 compatible (v1)" + '\n')
            f.write("fastqc version: " + os.popen("fastqc --version").read().strip() + '\n')
            f.write("bowtie2 version: " + os.popen("bowtie2 --version").read().strip() + '\n')
            f.write("samtools version: " + os.popen("samtools --version").read().strip() + '\n')
            f.write("picard version: 2.20.2-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
            f.write("bedtools version: " + os.popen("bedtools --version").read().strip() + '\n\n')
            f.write("PARAMETERS\n")
            f.write("##########\n")
            f.write("genome reference for aligning: " + config["genomeRef"] + '\n')
            f.write("blacklist for filtering: " + config["blacklist"] + '\n')
            f.write("map quality threshold for filtering: " + config["mapq"] + '\n')
            f.write("align command: " + align + '\n\n')
            f.write("SUMMARY\n")
            f.write("#######\n")
            # summary stats over the samples
            with open(params.temp, "w") as g:
                g.write("SAMPLE\tRAW_READ_PAIRS\tPERCENT_ALIGNED\tESTIMATED_LIBRARY_SIZE\tPERCENT_DUPLICATED\tPERCENT_MITOCHONDRIAL\tREAD_PAIRS_POST_FILTER\tPEAK_INSERTIONS_TSS\tMAX_MYCOPLASMA_MAP\n")
                for ftp in params.files:
                    g.write(ftp + '\t')
                    g.write(os.popen("awk '{{if (FNR == 1) print $1}}' " + ftp + "/adapter_trim.log").read().strip() + '\t')
                    g.write(os.popen("awk '{{if (FNR == 15) print $1}}' " + ftp + "/*.alignlog").read().strip() + '\t')
                    g.write(os.popen("""awk '{{if (FNR == 8) print $11}}' """ + ftp + "/dups.log").read().strip() +'\t')
                    g.write(os.popen("""awk '{{if (FNR == 8) dec=$10}}END{{printf("%.2f%",100*dec)}}' """ + ftp + "/dups.log").read().strip() +'\t')
                    shell("samtools idxstats " + ftp + "/*trim.st.bam > " + ftp + "/" + ftp + ".idxstats.dat")
                    g.write(os.popen("""awk '{{sum+= $3; if ($1 == "chrM") mito=$3}}END{{printf("%.2f%",100*mito/sum) }}' """ + ftp + "/" + ftp + ".idxstats.dat").read().strip() +'\t')
                    g.write(os.popen("samtools idxstats " + ftp + """/*.st.all*rmdup.bam | awk '{{s+=$3}} END{{printf("%i", s/2)}}'""").read().strip() +'\t')
                    if os.path.exists(config["TSS"]):
                        g.write(os.popen("sort -nrk1,1 " + ftp + """/*RefSeqTSS | head -1 | awk '{{printf("%.3f", $1)}}' """).read().strip() +'\t')
                    else:
                        g.write("NA" +'\t')
                    g.write(os.popen("""awk 'index($1, "Mycoplasma")' """ + ftp + "/*R1*trim_screen.txt " + """| awk '{{printf("%.2f%\\n", 100*($2-$3)/$2)}}' """ + "| sort -nrk1,1 | head -1").read().strip() + '\n')
                    # finish clean up by moving index file
                    shell("mv " + ftp + "/*.st.bam.bai " + ftp + "/00_source/") 
                    shell("mv " + ftp + "/" + ftp + ".idxstats.dat " + ftp + "/00_source/")
        # append summary log to rest of summary
        shell("cat {params.temp} | column -t >> {output}")
        shell("rm {params.temp}")


