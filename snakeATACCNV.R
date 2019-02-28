## Identification of CNV regions from background ATACseq data

## Load packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Rsubread", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rsubread)
library(genomation)

## Generate hg38 Granges
nContigs = length(seqnames(Hsapiens))
hg38 = GRanges(seqnames=seqnames(Hsapiens),
             ranges=IRanges(start=rep(1, nContigs), end=seqlengths(Hsapiens)),
             strand=rep("*", nContigs))
##
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")


## Generate windows
## This will generate a GRangesList, where each element is a GRanges object with the windows for one chromosome
## Individual GRanges can be extracted with "chr1 <- hg38Windows[[1]]" etc.
#hg38Tiles <- tile(gr, width = 1000000)
hg38Windows <- slidingWindows(hg38, width=1000000L, step=500000L)


## Calculate the total number of reads in each window



## Load the called peaks as a GRanges object
bedFile <- "C:\\Users\\Jordan\\Documents\\atac\\SNU61-WT-01.all_peaks.narrowPeak"
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Remove the entries that are not on standard chromosomes (chr1-22, X, Y), often helps prevent problems downstream
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")


## Calculate the total number of reads in each peak


## For each window, find which peaks overlap with it
peakGeneOverlaps <- findOverlaps(extPeaks, hg38Genes)
