#! usr/bin/env python 
# npagane | 190701 | risca lab | helper file for rules

import os
import numpy as np
import datetime
import sys
import subprocess

################################
# parameters and functions
################################

# this is for ease of development
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/fastq2bam"
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/fastq2bam" 

# defaults for parameters set in fastq2bam.py exectuable file

# determine sample names and sample numbers from the working directory
def findFiles(fastqDir, samp): 
    WHOLEFILES = []
    for base, dirs, files in os.walk(fastqDir + "/"):
        for fastq in files:
            if fastq.endswith(".fastq.gz") and samp == fastq.split("_")[0]:
                tmp = fastq.split(".fastq.gz")[0]
                WHOLEFILES.append(tmp.split('_'))
    return WHOLEFILES

# generate structure of expected files 
def customFileExpand(ext, fastqDir, samp, dir = ''):
    WHOLEFILES = findFiles(fastqDir, samp)
    strout = []
    parts = set(WHOLEFILES[0])
    for i in range(1, len(WHOLEFILES)):
        parts = set(parts & set(WHOLEFILES[i]))
    parts=list(parts)
    uniq, inds = np.unique(WHOLEFILES[0], return_index = True)
    if dir == '':
        ftp = samp + '/'
    else:
        ftp = samp + '/' + dir + '/'
    for i in range(len(inds)):
        ind = inds.tolist().index(i)
        if uniq[ind] in parts:
            ftp += uniq[ind] + '_'
    ftp = ftp[0:-1] # remove last character to then add extension
    ftp += ext
    strout.append(ftp)
    return strout

# conditional expand function upon 2 conditions for inputs/output
def conditionalExpand_2(condition1, condition2, truetrue, truefalse, falsetrue, falsefalse):
    if condition1:
        if condition2:
            fpt = truetrue
        else:
            fpt = truefalse
    else:
        if condition2:
            fpt = falsetrue
        else:
            fpt = falsefalse
    return fpt

# determine lanes that a certain sample was run through 
def determine_lanes(fastqDir, samp):
    lanes = []
    for reads in findFiles(fastqDir, samp):
        for part in reads:
            if 'L00' in part:
                lanes.append('_L00' + part.split('L00')[1][0] + '_')
    if len(lanes) == 0:
            lanes.append('_')
    return list(set(lanes))

################################
# combine insert plots and summary (7)
################################
def fastq2bamSummary(sampleTxt, genomeRef, blacklist, mapq, TSS, fastqDir):
    output="fastq2bamRunSummary.log"
    align = "(bowtie2 -X2000 -p{threads} -x {config[genomeRef]} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}"
    files = []
    with open(sampleTxt, 'r') as ftp:
        for line in ftp:
            files.append(line.strip())
    temp = "tempSummary_fastq.log"
    # make nice pdf of insert distributions
    os.system("Rscript " + exeDir + "/scripts/plotisds_v2.R " + sampleTxt + " hist_data_withoutdups")
    print('\n###########################')
    print('fastq2bam pipeline complete')
    print('\n###########################')
    # get git commit of current pipeline
    pcommit = subprocess.Popen(["git rev-parse", "HEAD"], cwd="/rugpfs/fs0/risc_lab/store/risc_soft/fastq2bam")
    pcommit.wait()
    with open(output, "w") as f:
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
        f.write("SOFTWARE\n")
        f.write("########\n")
        f.write("pipeline version (fastq2bam git commit): " + pcommit)
        f.write("python version: " + str(sys.version_info[0]) + '\n')
        f.write("pyadapter_trim version: python3 compatible (v1)" + '\n')
        f.write("fastqc version: " + os.popen("fastqc --version").read().strip() + '\n')
        f.write("bowtie2 version: " + os.popen("bowtie2 --version").read().strip() + '\n')
        f.write("samtools version: " + os.popen("samtools --version").read().strip() + '\n')
        f.write("picard version: 2.20.2-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
        f.write("bedtools version: " + os.popen("bedtools --version").read().strip() + '\n\n')
        f.write("PARAMETERS\n")
        f.write("##########\n")
        f.write("genome reference for aligning: " + genomeRef + '\n')
        f.write("blacklist for filtering: " + blacklist + '\n')
        f.write("map quality threshold for filtering: " + mapq + '\n')
        f.write("align command: " + align + '\n\n')
        f.write("SUMMARY\n")
        f.write("#######\n")
        # summary stats over the samples
        with open(temp, "w") as g:
            g.write("SAMPLE\tRAW_READ_PAIRS\tPERCENT_ALIGNED\tESTIMATED_LIBRARY_SIZE\tPERCENT_DUPLICATED\tPERCENT_MITOCHONDRIAL\tREAD_PAIRS_POST_FILTER\tPEAK_INSERTIONS_TSS\tMAX_MYCOPLASMA_MAP\n")
            for ftp in files:
                g.write(ftp + '\t')
                raw_read=0
                for lane in determine_lanes(fastqDir, ftp):
                    raw_read+=int(os.popen("awk '{{if (FNR == 1) print $1}}' " + ftp + "/*" + lane + "*adapter_trim.log").read().strip())
                g.write(str(raw_read) + '\t')
                palign=0
                for lane in determine_lanes(fastqDir, ftp):
                    palign+=float(os.popen("awk '{{if (FNR == 15) print $1}}' " + ftp + "/*" + lane + "*.alignlog").read().strip()[0:-1])/len(determine_lanes(fastqDir, ftp))
                g.write(str(np.round(palign, 2)) + '%\t')
                g.write(os.popen("""awk '{{if (FNR == 8) print $11}}' """ + ftp + "/dups.log").read().strip() +'\t')
                g.write(os.popen("""awk '{{if (FNR == 8) dec=$10}}END{{printf("%.2f%",100*dec)}}' """ + ftp + "/dups.log").read().strip() +'\t')
                os.system("samtools idxstats " + ftp + "/*trim.st.bam > " + ftp + "/" + ftp + ".idxstats.dat")
                g.write(os.popen("""awk '{{sum+= $3; if ($1 == "chrM") mito=$3}}END{{printf("%.2f%",100*mito/sum) }}' """ + ftp + "/" + ftp + ".idxstats.dat").read().strip() +'\t')
                g.write(os.popen("samtools idxstats " + ftp + """/*.st.all*rmdup.bam | awk '{{s+=$3}} END{{printf("%i", s/2)}}'""").read().strip() +'\t')
                if os.path.exists(TSS):
                    g.write(os.popen("sort -nrk1,1 " + ftp + """/*RefSeqTSS | head -1 | awk '{{printf("%.3f", $1)}}' """).read().strip() +'\t')
                else:
                    g.write("NA" +'\t')
                mycoplasma=0
                for lane in determine_lanes(fastqDir, ftp):
                    mycoplasma+=float(os.popen("""awk 'index($1, "Mycoplasma")' """ + ftp + "/*" + lane + "R1*trim_screen.txt " + """| awk '{{printf("%.2f\\n", 100*($2-$3)/$2)}}' """ + "| sort -nrk1,1 | head -1").read().strip())/len(determine_lanes(fastqDir, ftp))
                g.write(str(np.round(mycoplasma, 2)) + '%\n')
                # finish clean up by moving index file
                os.system("mv " + ftp + "/*.st.bam.bai " + ftp + "/00_source/")
                os.system("mv " + ftp + "/" + ftp + ".idxstats.dat " + ftp + "/00_source/")
    # append summary log to rest of summary
    os.system("cat " + temp + " | column -t >> " + output)
    os.system("rm " + temp)
    return


################################
# success and summary (7)
################################
def ATACseqSummary(sampleTxt, chromSize):
    output="ATACseqRunSummary.log"
    callpeak = "macs2 callpeak -f BAM -t {input} -n {params} -B --SPMR --nomodel --shift -37 --extsize 73 --nolambda --keep-dup all --call-summits --slocal 10000" # or -75 150
    bam2bg = "bedtools genomecov -ibam {input} -5 -bg -g {config[chromSize]} > {output.bg}" # need to update this whenever the command in the fastq2bam.smk file changes, dont like this but better than makign the command blind (placing it here) and then importing it
    files = []
    with open(sampleTxt, 'r') as ftp:
        for line in ftp:
            files.append(line.strip())
    temp = "tempSummary_atac.log"
    print('\n###########################')
    print('ATAC-seq pipeline complete')
    print('\n###########################')
    # get git commit of current pipeline
    pcommit_fast = subprocess.Popen(["git rev-parse", "HEAD"], cwd="/rugpfs/fs0/risc_lab/store/risc_soft/fastq2bam")
    pcommit_fast.wait()
    pcommit_atac = subprocess.Popen(["git rev-parse", "HEAD"], cwd="/rugpfs/fs0/risc_lab/store/risc_soft/ATACseq")
    pcommit_atac.wait()
    with open(output, "w") as f:
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n\n')
        f.write("SOFTWARE\n")
        f.write("########\n")
        f.write("pipeline version (fastq2bam git commit): " + pcommit_fast)
        f.write("pipeline version (ATACseq git commit): " + pcommit_atac)
        f.write("python version: " + str(sys.version_info[0]) + '\n')
        f.write("bedtools version: " + os.popen("bedtools --version").read().strip() + '\n')
        f.write("macs2 version: 2.1.2 <in macs2_python2.yml conda env>\n") # must update if macs2_python2 conda env is updated
        f.write("ucsc tools version: 2 (conda 332)\n\n") # must update if new version ever downloaded (shouldnt bc software dependencies)
        f.write("PARAMETERS" + '\n')
        f.write("##########\n")
        f.write("chromosome sizes: " + chromSize + '\n')
        f.write("peak call command: " + callpeak + '\n')
        f.write("bam to bedgraph command: " + bam2bg + '\n\n')
        f.write("SUMMARY\n")
        f.write("#######\n")
        with open(temp, "w") as g:
            g.write("SAMPLE\tFRACTION_READS_IN_PEAKS\n")
            for ftp in files:
                g.write(ftp + '\t')
                g.write(os.popen("""calc(){ awk "BEGIN { print "$*" }"; }; """ +  "den=`samtools view -c " + ftp + "/*.rmdup.atac.bam`; num=`bedtools sort -i " +
                ftp + "/peakCalls/*_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a " +ftp + "/*.rmdup.atac.bam -b stdin -ubam "
                    + "| samtools view -c`; calc $num/$den " + """| awk '{{printf("%.2f", $1)}}' """).read().strip() + '\n')
    # append summary log to rest of summary
    os.system("cat " + temp + " | column -t >> " + output)
    os.system("rm " + temp)
    return

#########################
# summary stats
if __name__ == '__main__':
    run=int(sys.argv[1])
    slurm=sys.argv[2]
    sampleTxt=sys.argv[3]
    g=open(sampleTxt, 'r')
    numSamples=len(g.readlines())
    g.close()
    sums=0
    for i in range(1,numSamples+1):
        f=open('slurm-' + slurm + '_' + str(i) + '.out', 'r')
        dat=f.readlines()
        if ('Nothing to be done.' in ' '.join(dat)):
            sums+=1
        f.close()
    if (sums < numSamples):
        if (not run):
            genomeRef=sys.argv[4]
            blacklist=sys.argv[5]
            mapq=sys.argv[6]
            TSS=sys.argv[7]
            fastqDir=sys.argv[8]
            fastq2bamSummary(sampleTxt, genomeRef, blacklist, mapq, TSS, fastqDir)
        else:
            chromSize=sys.argv[4]
            ATACseqSummary(sampleTxt, chromSize)
    else:
        print('data was not processed further')
