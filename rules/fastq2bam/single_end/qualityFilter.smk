#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# map quality filter (7)
################################

rule qualityFilter:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.qft.bam"
    run:
        print("Removing low quality reads")
        shell("samtools view -bh -q {config[mapq]} {input} -o {output}")
