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
exeDir="/rugpfs/fs0/risc_lab/store/risc_soft/pipeSeq"
#exeDir="/rugpfs/fs0/risc_lab/store/npagane/pipeSeq" 

# BOWTIE2 ALIGN COMMAND
align = "(bowtie2 -X2000 -p{threads} -x {config[genomeRef]} -1 {input.unzip1} -2 {input.unzip2} | samtools view -bS - -o {output.bam}) 2>{output.alignLog}"

# MACS2 PEAK CALL COMMAND
callpeak = "macs2 callpeak -f BAM -t {input} -n {params} -B --SPMR --nomodel --shift -37 --extsize 73 --nolambda --keep-dup all --call-summits --slocal 10000" # or -75 150

# BAMTOBEDGRAPH COMMAND
bam2bg = "bedtools genomecov -ibam {input} -5 -bg -g {config[chromSize]} > {output.bg}"

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
    # remove lane information if it remains
    newparts = []
    for part in parts:
        if 'L00' not in part:
            newparts.append(part)
    parts=newparts
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

# summary stats over the samples
def sampleSummaryStats(temp, files, fastqDir, TSS=None):
    with open(temp, "w") as g:
        if temp == "tempSummary_atac.log":
            g.write("SAMPLE\tRAW_READ_PAIRS\tPERCENT_ALIGNED\tESTIMATED_LIBRARY_SIZE\tPERCENT_DUPLICATED\tPERCENT_MITOCHONDRIAL\tREAD_PAIRS_POST_FILTER\tAVG_MYCOPLASMA_MAP\tADAPTER_MAP\tPEAK_INSERTIONS_TSS\tFRACTION_READS_IN_PEAKS\n")
        else:
            g.write("SAMPLE\tRAW_READ_PAIRS\tPERCENT_ALIGNED\tESTIMATED_LIBRARY_SIZE\tPERCENT_DUPLICATED\tPERCENT_MITOCHONDRIAL\tREAD_PAIRS_POST_FILTER\tAVG_MYCOPLASMA_MAP\tADAPTER_MAP\n")
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
            mycoplasma=0
            for lane in determine_lanes(fastqDir, ftp): # assume r1 is same as r2
                mycoplasma+=float(os.popen("""awk 'index($1, "plasma")' """ + ftp + "/*" + lane + "R1*trim_screen.txt " + """| awk '{{printf("%.2f\\n", 100*($2-$3)/$2)}}' """
                    + "| sort -nrk1,1 | head -1").read().strip())/len(determine_lanes(fastqDir, ftp))
            g.write(str(np.round(mycoplasma, 2)) + '%\t')
            adapter=0
            for lane in determine_lanes(fastqDir, ftp): # assume r1 is same as r2
                adapter+=float(os.popen("""awk 'index($1, "Adapters")' """ + ftp + "/*" + lane + "R1*trim_screen.txt " + """| awk '{{printf("%.2f\\n", 100*($2-$3)/$2)}}' """ 
                    + "| sort -nrk1,1 | head -1").read().strip())/len(determine_lanes(fastqDir, ftp))
            g.write(str(np.round(adapter, 2)) + '%')
            if temp == "tempSummary_atac.log":
                g.write('\t')
                g.write(os.popen("sort -nrk1,1 " + ftp + """/*RefSeqTSS.log | head -1 | awk '{{printf("%.3f", $1)}}' """).read().strip() +'\t')
                g.write(os.popen("""calc(){ awk "BEGIN { print "$*" }"; }; """ +  "den=`samtools view -c " + ftp + "/*.rmdup.atac.bam`; num=`bedtools sort -i " 
                    + ftp + "/peakCalls/*_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a " +ftp + "/*.rmdup.atac.bam -b stdin -ubam "
                    + "| samtools view -c`; calc $num/$den " + """| awk '{{printf("%.2f", $1)}}' """).read().strip())
                # finish clean up into intermediates directory
                os.system("if [ ! -d " + ftp + "/intermediates ]; then mkdir " + ftp + "/intermediates; fi")
                os.system("mv " + ftp + "/*.atac.bam " + ftp + "/intermediates/")
                os.system("rm -r " + ftp + "/.conda") # remove the conda software directory (NECESSARY TO SAVE DISK QUOTA SPACE)
            g.write('\n')
            # finish clean up into intermediates directory
            os.system("if [ ! -d " + ftp + "/intermediates ]; then mkdir " + ftp + "/intermediates; fi")
            os.system("mv " + ftp + "/*.trim.fastq.gz " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.trim.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.st.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.all.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.chrM.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.blft.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.qft.bam " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/*.bai " + ftp + "/intermediates/")
            os.system("mv " + ftp + "/" + ftp + ".idxstats.dat " + ftp + "/intermediates/")
            # remove the .snakemake folder upon successful completion
            os.system("rm -rf ./snakemake")
    return


################################
# combine insert plots and summary (7)
################################
def fastq2bamSummary(sampleTxt, invocationCommand, fastqDir, genomeRef, blacklist, mapq):
    output="fastq2bamRunSummary.log"
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
    pcommit = subprocess.Popen(["git", "rev-parse", "HEAD"], cwd="/rugpfs/fs0/risc_lab/store/risc_soft/pipeSeq", stdout=subprocess.PIPE)
    pcommit.wait()
    with open(output, "w") as f:
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n')
        f.write("invocation command: " + invocationCommand + '\n\n')
        f.write("SOFTWARE\n")
        f.write("########\n")
        f.write("pipeline version (pipeSeq git commit): " + pcommit.stdout.read().decode().strip() + '\n')
        f.write("python version: " + '.'.join(np.asarray(sys.version_info, dtype='str')[0:3]) + '\n')
        f.write("pyadapter_trim version: python3 compatible (v1 or v2 (same but 4x faster))" + '\n')
        f.write("fastqc version: " + os.popen("fastqc --version").read().strip() + '\n')
        f.write("bowtie2 version: " + os.popen("bowtie2 --version").read().strip() + '\n')
        f.write("samtools version: " + os.popen("samtools --version").read().strip() + '\n')
        f.write("picard version: 2.20.3-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
        f.write("bedtools version: " + os.popen("bedtools --version").read().strip() + '\n\n')
        f.write("PARAMETERS\n")
        f.write("##########\n")
        f.write("genome reference for aligning: " + genomeRef + '\n')
        f.write("blacklist for filtering: " + blacklist + " (file exists: " + str(os.path.exists(blacklist)) + ")" +'\n')
        f.write("map quality threshold for filtering: " + mapq + '\n')
        f.write("align command: " + align + '\n\n')
        f.write("SUMMARY\n")
        f.write("#######\n")
        # summary stats over the samples
        sampleSummaryStats(temp, files, fastqDir)
   # append summary log to rest of summary
    os.system("cat " + temp + " | column -t >> " + output)
    os.system("rm " + temp)
    return


################################
# success and summary (7)
################################
def ATACseqSummary(sampleTxt, invocationCommand, fastqDir, genomeRef, blacklist, mapq, TSS, chromSize):
    output="ATACseqRunSummary.log"
    files = []
    with open(sampleTxt, 'r') as ftp:
        for line in ftp:
            files.append(line.strip())
    temp = "tempSummary_atac.log"
    # make nice pdf of insert distributions
    os.system("Rscript " + exeDir + "/scripts/plotisds_v2.R " + sampleTxt + " hist_data_withoutdups")
    print('\n###########################')
    print('ATAC-seq pipeline complete')
    print('\n###########################')
    # get git commit of current pipeline
    pcommit = subprocess.Popen(["git", "rev-parse", "HEAD"], cwd="/rugpfs/fs0/risc_lab/store/risc_soft/pipeSeq", stdout=subprocess.PIPE)
    pcommit.wait()
    with open(output, "w") as f:
        f.write('user: ' + os.environ.get('USER') + '\n')
        f.write('date: ' + datetime.datetime.now().isoformat() + '\n')
        f.write("invocation command: " + invocationCommand + '\n\n')
        f.write("SOFTWARE\n")
        f.write("########\n")
        f.write("pipeline version (pipeSeq git commit): " + pcommit.stdout.read().decode().strip() + '\n')
        f.write("python version: " + '.'.join(np.asarray(sys.version_info, dtype='str')[0:3])  + '\n')
        f.write("pyadapter_trim version: python3 compatible (v1 or v2 (same but 4x faster))" + '\n')
        f.write("fastqc version: " + os.popen("fastqc --version").read().strip() + '\n')
        f.write("bowtie2 version: " + os.popen("bowtie2 --version").read().strip() + '\n')
        f.write("samtools version: " + os.popen("samtools --version").read().strip() + '\n')
        f.write("picard version: 2.20.3-SNAPSHOT" + '\n') # DONT LIKE THIS but the following wont work #+ os.popen("picard SortSam --version").read() + '\n')
        f.write("bedtools version: " + os.popen("bedtools --version").read().strip() + '\n')
        f.write("macs2 version: 2.1.2 <in macs2_python2.yaml conda env>\n") # must update if macs2_python2 conda env is updated
        f.write("ucsc tools version: 2 (conda 332)\n\n") # must update if new version ever downloaded (shouldnt bc software dependencies)
        f.write("PARAMETERS" + '\n')
        f.write("##########\n")
        f.write("genome reference for aligning: " + genomeRef + '\n')
        f.write("blacklist for filtering: " + blacklist + " (file exists: " + str(os.path.exists(blacklist)) + ")" +'\n')
        f.write("map quality threshold for filtering: " + mapq + '\n')
        f.write("TSS bed file for insertion calculation: " + TSS + '\n')
        f.write("chromosome sizes: " + chromSize + '\n')
        f.write("align command: " + align + '\n')
        f.write("peak call command: " + callpeak + '\n')
        f.write("bam to bedgraph command: " + bam2bg + '\n\n')
        f.write("SUMMARY\n")
        f.write("#######\n")
        # summary stats over the samples
        sampleSummaryStats(temp, files, fastqDir, TSS)
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
    invocationCommand=sys.argv[4]
    fastqDir=sys.argv[5]
    genomeRef=sys.argv[6]
    blacklist=sys.argv[7]
    mapq=sys.argv[8]
    sums=0
    for i in range(1,numSamples+1):
        f=open('slurm-' + slurm + '_' + str(i) + '.out', 'r')
        dat=f.readlines()
        if ('Nothing to be done.' in ' '.join(dat)):
            sums+=1
        f.close()
    if (sums < numSamples):
        if (not run):
            fastq2bamSummary(sampleTxt, invocationCommand, fastqDir, genomeRef, blacklist, mapq)
        else:
            TSS=sys.argv[9]
            chromSize=sys.argv[10]
            ATACseqSummary(sampleTxt, invocationCommand, fastqDir, genomeRef, blacklist, mapq, TSS, chromSize)
    else:
        print('data was not processed further')
