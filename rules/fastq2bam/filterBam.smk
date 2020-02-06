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
        shell("""
              if [ `echo {config[genomeRef]} | grep hg'[0-9]\{{2\}}' | wc -l` != 0 ]
              then
                  num=22
                  mito="MT"
              elif [ `echo {config[genomeRef]} | grep mm'[0-9]\{{2\}}' | wc -l` != 0 ]
              then 
                  num=19
                  mito="M"
              else
                  echo "cannot determine original alignment genome for further filtering"
              fi
              if [ `samtools view -H {input} | grep chr | wc -l` == 0 ]
              then
                  samtools view -o {output} {input} `seq 1 $num | sed -e '\$aX'`
                  echo `seq 1 $num | sed -e '\$aX'`
                  samtools view -b {input} $mito > {params.chrM}
              else
                  samtools view -o {output} {input} `seq 1 $num | sed 's/^/chr/' | sed -e '\$achrX'`
                  echo `seq 1 $num | sed 's/^/chr/' | sed -e '\$achrX'`
                  samtools view -b {input} chr$mito > {params.chrM}
              fi
              """)
