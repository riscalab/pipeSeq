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
        chrM = config['sample'] + "/{pre_tag}_{post_tag}.trim.st.chrM.bam"
    run:
        shell("samtools index {input}")
        print("Removing reads from unwanted chromosomes and scaffolds")
        shell('chrs="";for i in {{1..100}}; do chrs+=" chr$i"; done; chrs+=" chrX"; samtools view -b {input} `echo $chrs` > {output}')
        shell("samtools view -b {input} chrM > {params.chrM}")
