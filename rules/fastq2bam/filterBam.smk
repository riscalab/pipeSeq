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
              # continue to the next section of the genome numbering is compatible with a simple numeric system
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
              elif [ `echo {config[genomeRef]} | grep GRCg7b | wc -l` != 0 ] || [ `echo {config[genomeRef]} | grep GRCg7w | wc -l` != 0 ]
              then
                  continue=false
                  mito="M"
		  samtools view -o {output} {input} `printf "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchr23\nchr24\nchr25\nchr26\nchr27\nchr28\nchr29\nchr30\nchr31\nchr32\nchr33\nchr34\nchr35\nchr36\nchr37\nchr38\nchr39\nchrW\nchrZ"`
                  echo `printf "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchr23\nchr24\nchr25\nchr26\nchr27\nchr28\nchr29\nchr30\nchr31\nchr32\nchr33\nchr34\nchr35\nchr36\nchr37\nchr38\nchr39\nchrW\nchrZ"`
                  samtools view -b {input} $mito > {params.chrM}
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
              elif [ `echo {config[genomeRef]} | grep HBV | wc -l` != 0 ]
              then 
                  continue=false
                  cp {input} {output}
                  echo "NO FILTERING AT THIS STEP. ONLY ONE CONTINUOUS DNA SEGMENT"
                  touch {params.chrM}
	      elif [ `echo {config[genomeRef]} | grep hfxDS2 | wc -l` != 0 ]
              then
                  continue=false
                  samtools view -o {output} {input} `printf "chr1\nchrpHV1\nchrpHV2\nchrpHV3\nchrpHV4"`
                  echo `printf "chr1\nchrpHV1\nchrpHV2\nchrpHV3\nchrpHV4"`
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
