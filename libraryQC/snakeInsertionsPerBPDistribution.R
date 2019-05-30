#### Code testing ####
bamPath <- "C:\\Users\\jsk33\\Desktop\\LNCaP-WT-01-REP1of1.u.bam"

####
cat("Calculating per bp insertion distribution for:", sampleName, "on chromosome", currentChr, "\n")

#### Set snakemake variables
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
dirPath <- snakemake@wildcards[["path"]]

####
cat("Loading libraries", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stats4))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))

#### Make the hg38 GRanges reference
cat("Generating hg38 reference genome", "\n")
hg38 <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
## subset to standard chromosomes only
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")
hg38 <- trim(hg38)

####
cat("Loading relevant reads", "\n")
bamFile <- BamFile(bamPath)
param <- ScanBamParam(which = hg38)
bamIn <- readGAlignments(bamFile, param = param)
##
cat("Converting reads to insertions", "\n")
grIn <- granges(bamIn)
grIn <- keepStandardChromosomes(grIn, pruning.mode="coarse")
grIn <- trim(grIn)
grIn <- resize(grIn, width = 1)
##
cat("Shifting insertions +4/-5 bp", "\n")
plusIdx <- which(strand(grIn) == "+")
minusIdx <- which(strand(grIn) == "-")
grPlus <- grIn[plusIdx]
grMinus <- grIn[minusIdx]
## Shift
grPlusShifted <- shift(grPlus, shift=4L)
grMinusShifted <- shift(grMinus, shift=-5L)
## Merge
grMerged <- c(grPlusShifted, grMinusShifted)
shiftedInsertions <- grMerged
## 
cat("Generating insertion matrix", "\n")
insRLE <- coverage(grMerged)
## Create a views object for the Rle list using the Granges sites data
insViews <- Views(insRLE, grIn)
## Convert to a vector
insMatrix <- as.matrix(insViews)

####
uniqueIns <- unique(insMatrix)
numUniqueIns <- length(uniqueIns)
##
insData <- matrix(data = NA, ncol = 2, nrow = numUniqueIns)
colnames(insData) <- c("Number of Insertions", "Times Observed in Genome")

##
for (a in 1:numUniqueIns){
  insData[a,1] <- uniqueIns[a]
  numObs <- length(which(insMatrix == uniqueIns[a]))
  insData[a,2] <- numObs
}


## Save the data
save(sampleTotalReads, file = outPath)