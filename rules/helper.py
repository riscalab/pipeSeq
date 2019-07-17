#! usr/bin/env python 
# npagane | 190701 | risca lab | helper file for rules

import os
import numpy as np

################################
# parameters and functions
################################

# defaults for parameters set in fastq2bam.py exectuable file

# determine sample names and sample numbers from the working directory
def findFiles(exclusion): 
    WHOLEFILES = {}
    # check for exclusion file
    exclude = []
    if os.path.exists(exclusion):
        with open(exclusion, 'r') as ftp:
            for line in ftp:
                exclude.append(line.strip())
    for base, dirs, files in os.walk("."):
        for fastq in files:
            if fastq.endswith(".fastq.gz") and not fastq.startswith("Undetermined") and "trim" not in fastq:
                tmp1 = fastq.split(".fastq.gz")[0]
                tmp2 = tmp1.split('_')[0]
                # check for exclusion
                if tmp2 not in exclude:
                    if tmp2 not in WHOLEFILES.keys():
                        WHOLEFILES[tmp2] = []
                    WHOLEFILES[tmp2].append(tmp1.split('_'))
    return WHOLEFILES

# generate strucutre of expected files 
def customFileExpand(ext, exclusion):
    WHOLEFILES = findFiles(exclusion)
    strout = []
    for sample in WHOLEFILES.keys():
        parts = list(set(WHOLEFILES[sample][0]) & set(WHOLEFILES[sample][1])) # may be up to 4 but has to have at least 2
        uniq, inds = np.unique(WHOLEFILES[sample][0], return_index = True) 
        ftp = sample + '/'
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
