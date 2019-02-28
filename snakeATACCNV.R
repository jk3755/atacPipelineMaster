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
library(Rsamtools)

## Generate hg38 Granges
nContigs = length(seqnames(Hsapiens))
hg38 = GRanges(seqnames=seqnames(Hsapiens),
               ranges=IRanges(start=rep(1, nContigs), end=seqlengths(Hsapiens)))
##
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")

## Generate windows
## This will generate a GRangesList, where each element is a GRanges object with the windows for one chromosome
## Individual GRanges can be extracted with "chr1 <- hg38Windows[[1]]" etc.
#hg38Tiles <- tile(gr, width = 1000000)
hg38Windows <- slidingWindows(hg38, width = 1000000L, step = 500000L)
hg38Windows <- unlist(hg38Windows)

## Calculate the total number of reads in each window

### testing
bamFile <- "/home/ubuntu2/atac/ls1034/wt01/preprocessing/12all/LS1034-WT-01.all.bam"

## Rsubread required that you add a metadata column 'id' to the GRanges object
hg38Windows@elementMetadata$id <- c(1:length(hg38Windows@ranges@start))

## Create the annotation object
annotHg38Windows <- createAnnotationFile(hg38Windows)

# Calculate the number of reads
countHg38Windows <- featureCounts(files = bamFile, nthreads = 20, annot.ext = annotHg38Windows, isPairedEnd = TRUE)


## with Rsamtools

#
quickBamFlagSummary(bamFile)

params <- ScanBamParam(which = hg38Windows)
aln <- countBam(bamFile, param = params)




## Load the called peaks as a GRanges object
bedFile <- "C:\\Users\\Jordan\\Documents\\atac\\SNU61-WT-01.all_peaks.narrowPeak"
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Remove the entries that are not on standard chromosomes (chr1-22, X, Y), often helps prevent problems downstream
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")


## Calculate the total number of reads in each peak


## For each window, find which peaks overlap with it
peakGeneOverlaps <- findOverlaps(extPeaks, hg38Genes)
