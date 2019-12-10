#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

################################
# trim adapters (1)
################################

rule trimAdapters:
    input:
        r1 = config['fastqDir'] + "/{pre_tag}_R1_{post_tag}.fastq.gz",
        r2 = config['fastqDir'] + "/{pre_tag}_R2_{post_tag}.fastq.gz"
    output:
        config['sample'] + "/{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        config['sample'] + "/{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    params:
        r1 = "{pre_tag}_R1_{post_tag}.trim.fastq.gz",
        r2 = "{pre_tag}_R2_{post_tag}.trim.fastq.gz"
    run:
        shell("""
              f1=`zcat {input.r1} | awk '$1 ~ /^+/' | wc -l`;
              f2=`zcat {input.r2} | awk '$1 ~ /^+/' | wc -l`;
              l1=`zcat {input.r1} | wc -l`;
              l2=`zcat {input.r2} | wc -l`;
              l14=`expr $l1 '/' 4`;
              l24=`expr $l2 '/' 4`;
              if [ $f1 == $l14 ] && [ $f2 == $l24 ]; then echo 'fast trim'; python {workflow.basedir}/scripts/pyadapter_trimP3V2.py -a {input.r1} -b {input.r2} > {config[sample]}/{wildcards.pre_tag}_adapter_trim.log;
              else echo 'slow trim'; python {workflow.basedir}/scripts/pyadapter_trimP3.py -a {input.r1} -b {input.r2} > {config[sample]}/{wildcards.pre_tag}_adapter_trim.log; fi
            """)
        shell("mv {params.r1} {config[sample]}/")
        shell("mv {params.r2} {config[sample]}/")
