#! usr/bin/env snakemake
## npagane | risca lab | ATACseq pipeline rules

import sys
sys.path.append('{workflow.basedir}')
import rules.helper

################################
# commands with custom flags
################################      

callpeak = helper.callpeak_ATACseq # go to helper file to see / edit callpeak command

################################
# call peaks to make bed (2)
################################

rule callPeakSummits:
    input:
        config['sample'] + "/{pre_tag}_{post_tag}{ext}.rmdup.atac.bam"
    output:
        expand(config['sample'] + "/peakCalls/{{pre_tag}}_{{post_tag}}{{ext,.*}}.rmdup.atac{end}", end=["_summits.bed", "_peaks.xls", "_peaks.narrowPeak", "_treat_pileup.bdg", "_control_lambda.bdg"])
    params:
        config['sample'] + "/peakCalls/{pre_tag}_{post_tag}{ext,.*}.rmdup.atac"
    conda:
        "../../envs/macs2_python2.yaml" # path relative to current file, not working directory
    shell:
        callpeak # macs2 callpeak command defined above

