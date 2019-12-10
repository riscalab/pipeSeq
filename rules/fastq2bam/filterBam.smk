#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# bam filtering with samtools (5)
################################

rule filterBam:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.all.bam"
    params:
        chrM = config['sample'] + "/{pre_tag}_{post_tag}.trim.st.chrM.bam",
        filterLog = "filtering.log"
    run:
        shell("samtools index {input}")
        shell("echo 'Removing reads from unwanted chromosomes and scaffolds'")
        shell('chrs="";for i in {{1..100}}; do chrs+=" chr$i"; done; chrs+=" chrX"; samtools view -b {input} `echo $chrs` > {output}')
        shell("samtools view -b {input} chrM > {params.chrM}")
        shell("echo 'Filtering file {input} by rules filterBam and filter_removeDups_and_enrichTSS' >> {config[sample]}_{params.filterLog}")
        shell("mv {config[sample]}_{params.filterLog} {config[sample]}/{params.filterLog}")
