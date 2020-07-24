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
              continue=true
              if [ `echo {config[genomeRef]} | grep hg38 | wc -l` != 0 ] || [ `echo {config[genomeRef]} | grep hg19 | wc -l` != 0 ]
              then
                  num=22
                  mito="M"
              elif [ `echo {config[genomeRef]} | grep mm10 | wc -l` != 0 ]
              then 
                  num=19
                  mito="MT"
              elif [ `echo {config[genomeRef]} | grep mm9 | wc -l` != 0 ]
              then 
                  num=19
                  mito="M"
              elif [ `echo {config[genomeRef]} | grep EF2 | wc -l` != 0 ]
              then
                  continue=false
                  mito="MT"
                  samtools view -o {output} {input} `printf "I\nII\nIII\n"`
                  echo `printf "I\nII\nIII\n"`
                  samtools view -b {input} $mito > {params.chrM}
              elif [ `echo {config[genomeRef]} | grep dm6 | wc -l` != 0 ]
              then 
                  continue=false
                  mito="M"
                  samtools view -o {output} {input} `printf "2L\n2R\n3L\n3R\n4\n"`
                  echo `printf "2L\n2R\n3L\n3R\n4\n"`
                  samtools view -b {input} $mito > {params.chrM}
              else
                  echo "cannot determine original alignment genome for further filtering"
              fi
              if [ `samtools view -H {input} | grep chr | wc -l` == 0 ] && [ "$continue" == true ]
              then
                  samtools view -o {output} {input} `seq 1 $num | sed -e '\$aX'`
                  echo `seq 1 $num | sed -e '\$aX'`
                  samtools view -b {input} $mito > {params.chrM}
              elif [ "$continue" == true ]
              then
                  samtools view -o {output} {input} `seq 1 $num | sed 's/^/chr/' | sed -e '\$achrX'`
                  echo `seq 1 $num | sed 's/^/chr/' | sed -e '\$achrX'`
                  samtools view -b {input} chr$mito > {params.chrM}
              fi
              """)
