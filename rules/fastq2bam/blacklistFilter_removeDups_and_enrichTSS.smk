#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

import os

################################
# blacklist and map quality filter, 
# remove duplicates, and TSS enrichment (6)
################################

rule filter_removeDups_and_enrichTSS:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.all.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.bam"
    params:
        filterLog = config['sample'] + "/filtering.log",
        histDupsLog = config['sample'] + "/hist_data_withdups.log",
        histDupsPDF = config['sample'] + "/hist_graphwithdups.pdf",
        dupsLog = config['sample'] + "/dups.log",
        histNoDupsLog = config['sample'] + "/hist_data_withoutdups.log",
        histNoDupsPDF = config['sample'] + "/hist_graphwithoutdups.pdf"
    run:
        # filter by blacklist if provided
        if os.path.exists(config["blacklist"]):
            shell("echo 'Removing blacklisted reads'")
            shell("bedtools intersect -v -abam {input} -b {config[blacklist]} -wa > {config[sample]}_temp.bam") # produces temp file
            ftp = "{config[sample]}/{wildcards.pre_tag}_{wildcards.post_tag}.trim.st.all.blft.bam"
            shell("samtools view -bh -f 0x2 {config[sample]}_temp.bam -o " + ftp)
            shell("rm {config[sample]}_temp.bam") # remove temp file
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
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histDupsLog} H={params.histDupsPDF} W=1000 STOP_AFTER=5000000")
        # remove duplicates
        shell("echo 'Removing duplicates'")
        tmp = ftp.split('.bam')[0] + '.rmdup.bam'
        shell("picard MarkDuplicates I=" + ftp + " O=" + tmp + " METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        ftp = tmp
        # histogram without duplicates
        shell("samtools index " + ftp)
        shell("echo 'Histogram without duplicates'")
        shell("picard CollectInsertSizeMetrics I=" + ftp + " O={params.histNoDupsLog} H={params.histNoDupsPDF} W=1000 STOP_AFTER=5000000")
        # TSS enrichment if provided
        if os.path.exists(config["TSS"]):
            tmp = ftp.split('.bam')[0] + '.RefSeqTSS'
            shell(workflow.basedir + "/scripts/pyMakeVplot_css_v01.py -a " + ftp + " -b {config[TSS]} -e 2000 -p ends -s 5 -v -u --atac -o " + tmp)
        else:
            shell("echo 'TSS BED file not provided. not creating TSS enrichment profile'")
        # cleanup directory
        shell("if [ ! -d {config[sample]}/intermediates ]; then mkdir {config[sample]}/intermediates; fi")
        shell("mv {config[sample]}/*.trim.fastq.gz {config[sample]}/intermediates/")
        shell("mv {config[sample]}/*.trim.bam {config[sample]}/intermediates/")
        shell("mv {config[sample]}/*.all.bam {config[sample]}/intermediates/")
        shell("mv {config[sample]}/*.chrM.bam {config[sample]}/intermediates/")
        shell("mv {config[sample]}/*.blft.bam {config[sample]}/intermediates/")
        shell("mv {config[sample]}/*.qft.bam {config[sample]}/intermediates/")
