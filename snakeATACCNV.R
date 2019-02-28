
## Load packages
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

## Generate hg38 Granges
nContigs = length(seqnames(Hsapiens))
hg38 = GRanges(seqnames=seqnames(Hsapiens),
             ranges=IRanges(start=rep(1, nContigs), end=seqlengths(Hsapiens)),
             strand=rep("*", nContigs))
##
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")


## Generate windows
#hg38Tiles <- tile(gr, width = 1000000)
hg38Windows <- slidingWindows(hg38, width=1000000L, step=500000L)

