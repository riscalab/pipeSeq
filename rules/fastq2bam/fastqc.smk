#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# QC of fastq files (2)
################################

rule fastqc:
    input:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    output:
        expand(config['sample'] + "/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", "_fastqc.zip", ".fastq"])
    params:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    threads: 4
    run:
        shell("fastqc -t {threads} -o {config[sample]} {input.r1} {input.r2}")
        shell("unpigz {params.r1} {params.r2}")
