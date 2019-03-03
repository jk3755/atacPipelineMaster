## Identification of CNV regions from background ATACseq data
## Load packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genomation)
library(Rsamtools)

##
bamFile <- "/home/ubuntu2/atac/ls1034/wt01/preprocessing/12all/LS1034-WT-01.all.bam"

## Generate hg38 Granges
nContigs = length(seqnames(Hsapiens))
hg38 = GRanges(seqnames=seqnames(Hsapiens),
               ranges=IRanges(start=rep(1, nContigs), end=seqlengths(Hsapiens)))
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")


## Generate windows
## This will generate a GRangesList, where each element is a GRanges object with the windows for one chromosome
## Individual GRanges can be extracted with "chr1 <- hg38Windows[[1]]" etc.
#hg38Tiles <- tile(gr, width = 1000000)
hg38Windows <- slidingWindows(hg38, width = 1000000L, step = 500000L)
hg38Windows <- unlist(hg38Windows)


##
params <- ScanBamParam(which = hg38Windows)
aln <- countBam(bamFile, param = params)


## Load the called peaks as a GRanges object
bedFile <- "/home/ubuntu2/atac/ls1034/wt01/preprocessing/13allpeaks/LS1034-WT-01.all_peaks.narrowPeak"
ls1034Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Remove the entries that are not on standard chromosomes (chr1-22, X, Y), often helps prevent problems downstream
ls1034Peaks <- keepStandardChromosomes(ls1034Peaks, pruning.mode="coarse")


## For each window, subtract number of reads in peaks from total to get background reads
windowBackgroundReads <- matrix(data = NA, nrow = 6165, ncol = 4)
colnames(windowBackgroundReads) <- c("window", "total reads", "peak reads", "background reads")

for (a in 1:240){
  
  cat(a, "\n")
  windowBackgroundReads[a,1] <- a
  windowBackgroundReads[a,2] <- aln[a,6]
  
  ## Get intersection of peaks with current window
  window <- hg38Windows[a]
  inter <- intersect(window, ls1034Peaks)
  
  if (length(inter) == 0){
    windowBackgroundReads[a,3] <- 0
    windowBackgroundReads[a,4] <- windowBackgroundReads[a,2]
    
  } else {
    
    ## Calc signal in the intersection
    sbp <- ScanBamParam(which = inter)
    peaksignal <- countBam(bamFile, param = sbp)
    
    ## 
    windowBackgroundReads[a,3] <- sum(peaksignal[,6])
    
    ##
    windowBackgroundReads[a,4] <- (windowBackgroundReads[a,2] - windowBackgroundReads[a,3])
  }
}


med <- median(windowBackgroundReads[,4])
me <- mean(windowBackgroundReads[,4])

