#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

################################
# align at insertion center (1a)
################################

rule ATACoffset:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.bam"
    params:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac.temp.bam"
    threads: 6
    run:
        shell("samtools index {input}") # suppress the pysam/htslib warning about the index file
        shell("alignmentSieve --numberOfProcessors {threads} --ATACshift --bam {input} -o {params}")
        shell("samtools sort -O bam -o {output} {params}") #sort (dont use picard it is too strict about bam formatting)
        shell("samtools index {output}") # regenerate index file
        shell("rm {params}")
