#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# sort bam files for filtering (4)
################################

rule sortBam:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.bam"
    params:
        "{pre_tag}_{post_tag}.trim.st.bam"
    run:
        shell("picard SortSam -Xmx8g I={input}  O={params}  SORT_ORDER=coordinate")
        shell("mv {params} {output}")
