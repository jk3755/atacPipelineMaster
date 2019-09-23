
## For generating GR intervals to feed into the generate seqbias model function

####
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genomation)

####
source("C:\\Users\\jsk33\\OneDrive\\Git\\atacPipelineMaster\\snakeResources\\scripts\\atacFunctions.R")

#### For reference genome hg38 ####

## Chr1
outputPath <- "C:\\Users\\jsk33\\OneDrive\\Git\\atacPipelineMaster\\snakeResources\\bed\\hg38.chr1.bed"
chrName <- "chr1"
#
hg38GR <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
com <- paste0("refGR <- hg38GR[seqnames(hg38GR) == '", chrName, "']")
eval(parse(text = com))
refGR[2] <- refGR
strand(refGR[1]) <- "+"
strand(refGR[2]) <- "-"
refGR@strand@values <- droplevels(refGR@strand@values)
refGR <- keepSeqlevels(refGR, chrName, pruning.mode = c("coarse"))
writeGRangesToBEDwithStrand(refGR, outputPath)
## Test BED import
testGR <- importBED(outputPath)


## Chr22
outputPath <- "C:\\Users\\jsk33\\OneDrive\\Git\\atacPipelineMaster\\snakeResources\\refgr\\hg38.chr22.GR.RData"
chrName <- "chr22"

hg38GR <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
com <- paste0("refGR <- hg38GR[seqnames(hg38GR) == '", chrName, "']")
eval(parse(text = com))
refGR[2] <- refGR
strand(refGR) <- "+"
strand(refGR[2]) <- "-"
save(hg38GR, file = outputPath)


#### For LNCAP custom genome assembly ####
