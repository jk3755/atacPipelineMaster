## Load libraries
install.packages("numbers")
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(MotifDb))
suppressMessages(library(numbers))

## Set snakemake variables
numMotifs <- 
  
## Load motifs
## Do this one by one for memory considerations
signalsMerged <- list()
totalSites <- 0
totalBp <- 0
motifWidths <- c()


##
load("C:/Users/jsk33/Desktop/SNU61-WT-01.CDX2.motif1.info.Rdata")
signals1 <- parsedSitesInfo[["bfPassPeakSignals"]]
sites1 <- parsedSitesInfo[["bfPassPeakSites"]]
numSites1 <- length(signals1[["+"]][,1])
motifWidth1 <- sites1@ranges@width[1]
motifWidths <- c(motifWidths, motifWidth1)
totalSites <- (totalSites + length(signals1[["+"]][,1]))

##
load("C:/Users/jsk33/Desktop/SNU61-WT-01.CDX2.motif2.info.Rdata")
signals2 <- parsedSitesInfo[["bfPassPeakSignals"]]
sites2 <- parsedSitesInfo[["bfPassPeakSites"]]
numSites2 <- length(signals2[["+"]][,1])
motifWidth2 <- sites2@ranges@width[1]
motifWidths <- c(motifWidths, motifWidth2)
totalSites <- (totalSites + length(signals2[["+"]][,1]))

## Must set the total bp to the length of the smallest vector
minWidth <- min(motifWidths)
totalBp <- min(motifWidths)+200

## Now, shave off values from vectors that are larger than the smallest motif
## Must retain positioning of the motif
mod <- mod((motifWidth2+1-minWidth),2)

## two potential patterns
## if mod is 0, can evenly split from either side
if (mod == 0)
  {
  ## regen
  }

## if mod is 1, need to remove 1 extra bp from one side (left side)
if (mod == 1){}


## To center the Granges, change the width of all sites to the size of smallest motif width
width(sites2) <- 8

## Once all Granges are centered together, can simply remove trailing values from the signals on right side
sig2 <- signals2[['+']]

##
mergePlus <- matrix(data = NA, nrow = totalSites, ncol = totalBp)
#
for (a in 1:numSites1){for (b in 1:totalBp){mergePlus[a,b] <- signals1[["+"]][a,b]}}
for (a in (1+numSites1):(numSites1+numSites2)){for (b in 1:totalBp){mergePlus[a,b] <- signals2[["+"]][(a-numSites1),b]}}

## use rbind to merge the matrices by row. Note, must handle the bp mismatch problem prior to this
signalsMerged$'+' <- mergePlus

one <- signals1[["+"]]
two <- signals2[["+"]]
new <- rbind(one,two)



