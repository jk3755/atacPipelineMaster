## Load libraries
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(MotifDb))

## Set snakemake variables
numMotifs <- 
  
## Load motifs
## Do this one by one for memory considerations
signalsMerged <- list()
totalSites <- 0
totalBp <- 0
totalBp <- (totalBp + length(signals1[["+"]][1,]))

##
signals1 <- parsedSitesInfo[["bfPassPeakSignals"]]
sites1 <- parsedSitesInfo[["bfPassPeakSites"]]
numSites1 <- length(signals1[["+"]][,1])
totalSites <- (totalSites + length(signals1[["+"]][,1]))

##
signals2 <- parsedSitesInfo[["bfPassPeakSignals"]]
sites2 <- parsedSitesInfo[["bfPassPeakSites"]]
numSites2 <- length(signals2[["+"]][,1])
totalSites <- (totalSites + length(signals2[["+"]][,1]))

##
mergePlus <- matrix(data = NA, nrow = totalSites, ncol = totalBp)
#
for (a in 1:numSites1){for (b in 1:totalBp){mergePlus[a,b] <- signals1[["+"]][a,b]}}
for (a in (1+numSites1):(numSites1+numSites2)){for (b in 1:totalBp){mergePlus[a,b] <- signals2[["+"]][(a-numSites1),b]}}

signalsMerged$'+' <- mergePlus

one <- signals1[["+"]]
two <- signals2[["+"]]
new <- rbind()



