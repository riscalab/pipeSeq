#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# The following program will trim PE reads

# Edited to add option for how many bases to align (leave it at n = 20, default)

##### IMPORT MODULES #####
# import necessary for python
import os
import re
import sys
import gzip
import Levenshtein
from Bio import SeqIO
from Bio import AlignIO
from optparse import OptionParser

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = str.maketrans('ATCGN', 'TAGCN') # NP EDIT change string to str (also not converting to bytes)
def reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-n", type=int,default=20,help="<n> accepts the number of bases to align to look for where to trim")
opts.add_option("-u", action="store_true", default=False, help="Print uncompressed output file")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# name input and outputs
p1_in = options.a
p2_in = options.b

# name outputs and print to working dir
p1_file = p1_in.split('/')[-1]
p2_file = p2_in.split('/')[-1]

#check for file type and open input file
append = p1_in.split('.')[-1]
if append == "fastq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = re.sub(".fastq", ".trim.fastq", p1_file)
    p2_out = re.sub(".fastq", ".trim.fastq", p2_file)
elif append == "fq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = re.sub(".fq", ".trim.fastq", p1_file)
    p2_out = re.sub(".fq", ".trim.fastq", p2_file)
elif append == "gz":
    p1_rds = gzip.open(p1_in,'rt') # NP EDIT from 'r' to 'rt' since 'r' defaults to 'rb'
    p2_rds = gzip.open(p2_in,'rt') 
    p1_out = re.sub(".fastq.gz", ".trim.fastq", p1_file)
    p2_out = re.sub(".fastq.gz", ".trim.fastq", p2_file)
else:
    sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0;
n=options.n  # match seq
mismatch=1  # only allow 0-1 mismatches for now

# initilize write files
if options.u == False:
    r1_write = gzip.open(p1_out+'.gz', 'wt') # NP EDIT to wt
    r2_write = gzip.open(p2_out+'.gz', 'wt') # same
elif options.u == True:
    r1_write = gzip.open(p1_out, 'wt') # NP EDIT to wt from w
    r2_write = gzip.open(p2_out, 'wt') # same

while 1:
    # read lines
    seqhead1 = p1_rds.readline()
    seqhead2 = p2_rds.readline()

    # break if at end of file
    if not seqhead1:
        break

    # load fastq into memory
    seq1 = p1_rds.readline().rstrip()
    seq2 = p2_rds.readline().rstrip()
    qualhead1 = p1_rds.readline()
    qualhead2 = p2_rds.readline()
    qual1 = p1_rds.readline().rstrip()
    qual2 = p2_rds.readline().rstrip()

    # align reads to themselves
    i = i+1  # total reads
    rc_seq2 = reverse_complement(seq2[0:n])
    idx = seq1.rfind(rc_seq2) # look for perfect match
    if idx > 0:
        j = j+1  # 0 mismatchs
    elif mismatch>0:
        hold = fuzz_align(rc_seq2,seq1,mismatch)  # else allow for mismatch
        if hold:
            idx,mis=hold
            if mis == 1:
                k=k+1  # 1 mismatch

    # trim reads if idx exist
    if idx > 0:
        # keep track on how much trimming
        tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
            
        # trim data
        seq1 = seq1[0:idx+n-1] # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
        seq2 = seq2[0:idx+n-1]
        qual1 = qual1[0:idx+n-1]
        qual2 = qual2[0:idx+n-1]
       
    # print read1
    r1_write.write(seqhead1)
    r1_write.write(seq1 + "\n")
    r1_write.write(qualhead1)
    r1_write.write(qual1 +"\n")

    # print read2
    r2_write.write(seqhead2)
    r2_write.write(seq2 + "\n")
    r2_write.write(qualhead2)
    r2_write.write(qual2 + "\n")


# close files to write the file
r1_write.close()
r2_write.close()
p1_rds.close()
p2_rds.close()

# give summary
print(str(i)+" sequences total") # change dumb python print statements NP 190624
print(str(j)+" sequences trimmed with 0 mismatches")
print(str(k)+" sequences trimmed with 1 mismatch")
# in case there is no trims NP EDIT (this was added to debug with smaller files that may not contain trims and because im scared of exception errors)
try:
   print(str(int(tot_b/(j+k)))+" mean number of bases trimmed for reads requiring trimming")
except ZeroDivisionError:
        print('divide by zero, there were no trimmed sequences')

