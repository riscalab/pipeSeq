#!/bin/bash
# npagane | risca lab | oct 2019 | CUTnTag pipeline wrapper to execute

##################
# SET PARAMETERS #
##################

# set the default arguments for optional parameters
genomeMap="hg38"
mapq=30
snakemake=""
exeDir=`readlink -f $0`
exeDir=`dirname $exeDir`

function print_usage {
    echo "USAGE: $0 [-c .] [-f /rugpfs/fs0/risc_lab/store/risc_data/runs/FASTQs] [-s sample.txt]"
    echo "  -c   PATH   path to where you want to process the sequencing data"
    echo "  -f   PATH   path to where your FASTQ files live"
    echo "OPTIONAL ARGUMENTS:"
    echo "  -g=hg38      STRING  genome to align FATSQ files to (i.e. hg38, hg19, mm9, mm10, etc.)"
    echo "  -b~hg38      PATH    path to where your blacklist file lives for filtering"
    echo "  -m=30        INT     map quality threshold for alignment"
    echo "  -p=''        STRING  additional snakemake flags and arguments for compilation"
    echo "  --singleend          if FASTQ files are single-end rather than paired-end sequenced"
    echo "  --help               print this usage text and exit"
}

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--singleend") set -- "$@" "-e foo";;
    "--help") set -- "$@" "-h bar" ;;
    *)        set -- "$@" "$arg"
  esac
done

# Parse short options
singleend="False"
OPTIND=1
while getopts c:f:g:b:m:p:e:h: option
do
case "${option}"
in
c) cwd=${OPTARG};; # working directory for analysis (REQUIRED, i.e. path/to/workingDirectory)
f) fastqDir=${OPTARG};; # directory with fastq files (REQUIRED, string, i.e. path/to/fastq)
g) genomeMap=${OPTARG};; # genome reference for alignment (OPTIONAL, string, OPTS: 'hg38', 'hg19', 'mm9', 'mm10')
b) blacklist=${OPTARG};; # blacklist for filtering (OPTIONAL, string, defaults to genomeMap blacklist OR overwrite with path/to/blacklist)
m) mapq=${OPTARG};; # the map quality threshold for alignment (OPTIONAL, int, i.e. 30)
p) snakemake=${OPTARG};; # any snakemake flags for compilation (OPTIONAL, string, i.e. "--snakemake unlock")
e) singleend="True";; # whether or not the sequenced data is single-end or not (i.e. paried-end)
h) print_usage; exit 0;;
esac
done

shift $(expr $OPTIND - 1) # remove options from positional parameters

# check for required arguments
if [ -z "$cwd" ] || [ -z "$fastqDir" ]
then
    echo "must specificy the desired working directory (-c) and fastq directory (-f)."
    echo ""
    print_usage
    exit
fi

# check genomeMap to align to correct genome with either assumed or specified blacklist
if [ "$genomeMap" != "hg38" ]
then 
    if [ "$genomeMap" == "hg19" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/blacklist/hg19-blacklist.v2.bed"
        fi
    elif  [ "$genomeMap" == "mm9" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm9/blacklist/mm9-blacklist.bed"
        fi
    elif [ "$genomeMap" == "mm10" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/blacklist/mm10-blacklist.v2.bed"
        fi
    elif [ "$genomeMap" == "dm6" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/dm6/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/dm6/blacklist/dm6-blacklist.v2.bed"
        fi
    elif [ "$genomeMap" == "EF2" ]
    then
        genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/S_pombe_EF2/genome/Sequence/Bowtie2Index/genome"
        if [ -z "$blacklist" ]
        then
            blacklist="None"
        fi
    else
        echo "unrecognized genome.\navailable genomes: hg38, hg19, mm10, mm9, dm6, EF2.\ntalk to nicole to get your genome on the cluster if not there.\n"
        exit
    fi
else
    genomeRef="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/Bowtie2Index/genome"
    if [ -z "$blacklist" ]
    then
        blacklist="/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/blacklist/hg38-blacklist.v2.bed"
    fi
fi
echo "alinging to $genomeRef with blacklist $blacklist"

# change working directory
cd $cwd

# set conda env
source activate fastq2bam

# make samples.txt file
echo "Generating samples.txt file" 
sampleText="samples.txt"
ls $fastqDir | awk -F'_S[-.0-9]*_' '{print $1}' | sort -u > $sampleText

# run CUTnTag 
numSamples=`wc -l $sampleText | awk '{print $1}' `
CUTnTag=$(sbatch -p hpc --array=1-$numSamples $exeDir/scripts/CUTnTag_snakemake.sh $cwd $fastqDir $sampleText $genomeRef $blacklist $mapq  $singleend $exeDir $snakemake)

# get job id
if ! echo ${CUTnTag} | grep -q "[1-9][0-9]*$"; then
   echo "Job(s) submission failed."
   echo ${CUTnTag}
   exit 1
else
   CUTnTagID=$(echo ${CUTnTag} | grep -oh "[1-9][0-9]*$")
fi

# summary stats for CUTnTag after successful completion
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
sbatch -p hpc --depend=afterok:$CUTnTagID --wrap="python $exeDir/rules/helper.py 2 $CUTnTagID $sampleText '$myInvocation' $fastqDir $genomeRef $blacklist $mapq $singleend $exeDir"
