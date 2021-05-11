#! usr/bin/env snakemake
## npagane | risca lab | ATAC pipeline rules

################################
# TSS enrichment (1b)
################################

rule enrichTSS:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.bam"
    output:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.RefSeqTSS.log"
    params:
        config['sample'] + "/{pre_tag}_{post_tag}{ext,.*}.RefSeqTSS"
    run:
        if os.path.exists(config['TSS']):
            shell(workflow.basedir + "/scripts/pyMakeVplot_css_v01.py -a {input} -b {config[TSS]} -e 2000 -p ends -s 6 -v -u --atac -o {params}")
        else:
            shell("touch {params}")
            shell("echo 0 > {output}")

