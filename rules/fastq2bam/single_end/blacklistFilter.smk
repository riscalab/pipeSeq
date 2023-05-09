#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# blacklist filter (6)
################################

rule blacklistFilter:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.st.all.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.blft.bam"
    run:
        print("Removing blacklisted reads")
        shell("bedtools intersect -v -abam {input} -b {config[blacklist]} -wa > {config[sample]}_temp.bam") # produces temp file
        shell("samtools view -bh {config[sample]}_temp.bam -o {output}")
        shell("rm {config[sample]}_temp.bam") # remove temp file
