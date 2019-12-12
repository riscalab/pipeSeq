#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

#########################################
# QC of fastq files and contamination (2)
#########################################

rule fastqcAndScreen:
    input:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    output:
        expand(config['sample'] + "/{{pre_tag}}_{run}_{{post_tag}}.trim{end}", run=["R1", "R2"], end=["_fastqc.html", "_fastqc.zip", "_screen.html", "_screen.png", "_screen.txt"]),
        unzip1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq",
        unzip2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq",
    params:
        r1 = config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    threads: 4
    run:
        shell("fastqc -t {threads} -o {config[sample]} {input.r1} {input.r2}")
        shell("unpigz {params.r1} {params.r2}")
        # go to /rugpfs/fs0/risc_lab/store/risc_soft/anaconda2/envs/fastq2bam/share/fastq-screen-0.13.0-0/fastq_screen.conf to make changes to fastq_screen
        shell("fastq_screen --aligner bowtie2 {output.unzip1} {output.unzip2}")
        shell("mv {wildcards.pre_tag}_R1_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
        shell("mv {wildcards.pre_tag}_R2_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
