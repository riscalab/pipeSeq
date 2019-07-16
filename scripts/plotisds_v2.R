#!/usr/bin/Rscript

# plotisds_v2.R
# Viviana Risca
# May 4, 2015
# Version 1.1 # June 1, 2015 changed hist_data_cleaned to hist_data_nodups
# Version 2.0 # June 25, moved libs.txt and hist_data_nodups_1linehdr into cmd line params
# edited by Nicole Pagane to remove extra lines at the beginning of hist log files

read.isd <- function(filename) {  
  # Read the ISD histogram file
  isdhist <- read.table(filename, header=TRUE, sep="\t", skip=10, col.names=c("fragment_length", "read_count"))
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
outfilename <- "isdplots"
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
  if ( max(isd$read_count) >= max_read) {
    max_read = max(isd$read_count)
  }
}

for (page in seq(from=1, to=numpages)) {
  
  pdf(sprintf(paste0(outfilename, "%03d", ".pdf"),page), paper="a4r", width=9, height=7)
  par(mfrow=calculate.plot.layout(plotsperpage))

  for (i in seq(from=minplot, to=maxplot)) {
    isdfilename <- make.hist.filename(libslist[i], logfileprefix, "log")
    isd <- read.isd(isdfilename)
    plot(isd$fragment_length, isd$read_count, type="l", xlab="length (bp)", ylab="reads", main=libslist[i], xlim=c(0,600), ylim=c(0, max_read))
  }
  
  minplot <- maxplot + 1
  maxplot <- min((maxplot + plotsperpage), numlibs)
  dev.off()
}


