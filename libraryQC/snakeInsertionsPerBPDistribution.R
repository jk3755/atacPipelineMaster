## Set snakemake variables
cat("Setting snakemake variables...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
dirPath <- snakemake@wildcards[["path"]]
currentChr <- snakemake@wildcards[["chr"]]

##
cat("Calculating per bp insertion distribution for:", sampleName, "on chromosome", currentChr, "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stats4))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(genomation))



#### Make the hg38 GRanges reference
hg38 <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
## subset to standard chromosomes only
hg38 <- keepStandardChromosomes(hg38, pruning.mode="coarse")
hg38 <- trim(hg38)

#### Subset hg38 Granges reference to current chromosome
hg38 <- hg38["chr1"]

##
bamFile <- BamFile(bamPath)
sampleTotalReads <- countBam(bamFile)
sampleTotalReads <- sampleTotalReads[[6]]
cat("Found", sampleTotalReads, "total reads in sample library", "\n")
      
## Save the data
save(sampleTotalReads, file = outPath)