#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# remove duplicates (8)
################################

rule removeDuplicates:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.bam"
    params:
        sorted = config['sample'] + "/{pre_tag}_{post_tag}{ext}.SORT.bam",
        histDupsLog = config['sample'] + "/hist_data_withdups.log",
        histDupsPDF = config['sample'] + "/hist_graphwithdups.pdf",
        dupsLog = config['sample'] + "/dups.log",
        histNoDupsLog = config['sample'] + "/hist_data_withoutdups.log",
        histNoDupsPDF = config['sample'] + "/hist_graphwithoutdups.pdf",
        stopAfter="5000000",
        width="1000"
    run:
        # sort one more time after all filtering
        shell("picard SortSam -Xmx8g I={input}  O={params.sorted}  SORT_ORDER=coordinate")
        shell("mv {params.sorted} {input}")
        # histogram with duplicates
        print("Histogram with duplicates")
        shell("picard CollectInsertSizeMetrics -Xmx8g I={input} O={params.histDupsLog} H={params.histDupsPDF} W={params.width} STOP_AFTER={params.stopAfter}")
        # remove duplicates
        print("Removing duplicates")
        shell("picard MarkDuplicates -Xmx8g I={input} O={output} METRICS_FILE={params.dupsLog} REMOVE_DUPLICATES=true")
        # histogram without duplicates
        shell("samtools index {output}")
        print("Histogram without duplicates")
        shell("picard CollectInsertSizeMetrics -Xmx8g I={output} O={params.histNoDupsLog} H={params.histNoDupsPDF} W={params.width} STOP_AFTER={params.stopAfter}")
