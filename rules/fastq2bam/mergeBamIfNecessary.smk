#! usr/bin/env snakemake
## npagane | risca lab | fastq2bam pipeline rules

import sys
sys.path.append('{workflow.basedir}')
import rules.helper

################################
# merge sample bams aross lanes (3.5)
################################

rule mergeBamIfNecessary:
    input:
        expand(config['sample'] + "/{{pre_tag}}{lanes}{{post_tag}}.trim.unmerged.bam", lanes=helper.determine_lanes(config["fastqDir"], config["sample"]))
    output:
        config['sample'] + "/{pre_tag}_{post_tag}.trim.bam"
    run:
        if (len(input) == 1):
            print("sample " + config['sample'] + " was sequenced in one lane")
            shell("mv {input} {output}")
        else:
            print('sample ' + config['sample'] + ' was sequenced in more than one lane; merging BAMs!')
            input_string = ''
            for inp in input:
                input_string += inp + ' '
            shell("samtools merge " + " {output} " + input_string)
