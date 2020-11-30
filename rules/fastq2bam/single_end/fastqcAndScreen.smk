#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

#########################################
# QC of fastq files and contamination (2)
#########################################

rule fastqcAndScreen:
    input:
        r = config['sample'] + "/{pre_tag}_{post_tag}.trim.fastq.gz"
    output:
        expand(config['sample'] + "/{{pre_tag}}_{{post_tag}}.trim{end}", end=["_fastqc.html", "_fastqc.zip", "_screen.html", "_screen.png", "_screen.txt"]),
        unzip = config['sample'] + "/{pre_tag}_{post_tag}.trim.fastq",
    params:
        r = config['sample'] + "/{pre_tag}_{post_tag}.trim.fastq.gz",
    threads: 4
    run:
        shell("fastqc -t {threads} -o {config[sample]} {input.r}")
        shell("unpigz {params.r}")
        # go to /rugpfs/fs0/risc_lab/store/risc_soft/anaconda2/envs/fastq2bam/share/fastq-screen-0.13.0-0/fastq_screen.conf to make changes to fastq_screen
        shell("fastq_screen --aligner bowtie2 {output.unzip}")
        shell("mv {wildcards.pre_tag}_{wildcards.post_tag}.trim_screen.* {config[sample]}/")
