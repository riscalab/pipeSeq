#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

################################
# create bigwig for tracks (5)
################################

rule bedGraph2bigWig:
    input:
        config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext}.rmdup.atac{ext2}.bdg"
    output:
        bw = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bw"
    params:
        qft = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bdg.clip",
        st = config['sample'] + "/{dir}/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac{ext2,.*}.bdg.st.clip"
    run:
        if (os.stat(input[0]).st_size > 50): 
            shell("bedtools slop -i {input} -g {config[chromSize]} -b 0 | bedClip stdin {config[chromSize]} {params.qft}")
            shell("LC_COLLATE=C sort -k1,1 -k2,2n {params.qft} > {params.st}")
            shell("bedGraphToBigWig {params.st} {config[chromSize]} {output.bw}")
            shell("rm {params.qft} {params.st}")
        else: 
            shell('cp {input} {output.bw}')
        # cleanup directory
        shell("if [ ! -d {config[sample]}/intermediates ]; then mkdir {config[sample]}/intermediates; fi")
        # this rule is duplicated for each sample, so actual clean up is in summary stats
