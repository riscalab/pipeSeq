#!/usr/bin/Rscript

# plotisds_v3.R
# Viviana Risca
# May 4, 2015
# Version 1.1 # June 1, 2015 changed hist_data_cleaned to hist_data_nodups
# Version 2.0 # June 25, moved libs.txt and hist_data_nodups_1linehdr into cmd line params
# edited by Nicole Pagane to ignore extra lines at the beginning of hist log files and to standardize y-axis
# Version 3: look at full histogram

read.isd <- function(filename) {  
  # Read the ISD histogram file
  isd <- read.table(filename)
  isd <- isd[isd>0]
  isd <- hist(isd, plot=FALSE, breaks=0:max(isd))
}

make.hist.filename <- function(lib, histprefix, ext) {
  temp_filename <- paste(lib,histprefix,sep="/") 
  filename <- paste(temp_filename, ext, sep=".")
  filename
}

calculate.plot.layout <- function(nplots) {
  ncol <- 4
  #  ncol <- ceiling(sqrt(nplots))
  nrow <- ceiling(nplots / ncol)
  layout <- c(nrow, ncol)
  layout
}

#### MAIN SCRIPT ####

# Read command line parameters
args<-commandArgs(TRUE)

if (length(args) == 0) {
        print("Usage: plotisds_v1.R libs.txt histogramlogfileprefix")
}

filename <- args[1]
logfileprefix <- args[2]

# Hard coded parameters
outfilename <- "isdplots_full_FLD"
plotsperpage <- 12  

# Read in ISD files from current directory and plot them onto the same PDF
libslist <- as.matrix(read.delim(filename, header=FALSE))
numlibs <- length(libslist)
numpages <- ceiling(numlibs/plotsperpage)

minplot <- 1
maxplot <- min(plotsperpage, numlibs)

# determine maximum read (y value) (NP EDIT)
max_read <- 1
for (i in seq(from=minplot, to=numlibs)) {
  isdfilename <- make.hist.filename(libslist[i], logfileprefix, "log")
  isd <- read.isd(isdfilename)
  if ( max(isd$counts) >= max_read) {
    max_read = max(isd$counts)
  }
}

for (page in seq(from=1, to=numpages)) {
  
  pdf(sprintf(paste0(outfilename, "%03d", ".pdf"),page), paper="a4r", width=9, height=7)
  par(mfrow=calculate.plot.layout(plotsperpage))

  for (i in seq(from=minplot, to=maxplot)) {
    isdfilename <- make.hist.filename(libslist[i], logfileprefix, "log")
    isd <- read.isd(isdfilename)
    plot(isd$mids, isd$counts, type="l", xlab="length (bp)", ylab="reads", main=libslist[i], xlim=c(0,600), ylim=c(0, max_read))
  }
  
  minplot <- maxplot + 1
  maxplot <- min((maxplot + plotsperpage), numlibs)
  dev.off()
}


