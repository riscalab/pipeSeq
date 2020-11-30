#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        r = config['fastqDir'] + "/{pre_tag}_{post_tag}.fastq.gz",
    output:
        r = config['sample'] + "/{pre_tag}_{post_tag}.trim.fastq.gz"
    params:
        r = "{pre_tag}_{post_tag}_trimmed.fq.gz",
        summary = "{pre_tag}_{post_tag}.fastq.gz_trimming_report.txt"
    run:
        shell("trim_galore {input.r}")
        shell("mv {params.r} {output.r}")
        shell("mv {params.summary} {config[sample]}/{wildcards.pre_tag}_adapter_trim.log")
