## See below link for reference information related to ChIPseeker package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("ChIPseeker", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("clusterProfiler", suppressUpdates = TRUE)
#biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene", suppressUpdates = TRUE)
#biocLite("org.Hs.eg.db", suppressUpdates = TRUE)
#biocLite("ReactomePA", suppressUpdates = TRUE)

##
cat("Loading libraries...", "\n")
library(ChIPseeker)
library(genomation)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##
cat("Setting snakemake vars...", "\n")
bedFile <- snakemake@input[[1]]
outputPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

## Data import
cat("Input peak file path:", bedFile, "\n")
narrowPeaks <- ChIPseeker::readPeakFile(bedFile)
## Subset to standard xsomes
narrowPeaks <- keepStandardChromosomes(narrowPeaks, pruning.mode="coarse")


## Coverage plots make plot of genome wide peak coverage
covplotPath <- gsub("annotations.done.txt", "repmerged.peakgenomecov.svg", outputPath)
covplotPath <- gsub("operations", "metrics", covplotPath)
cat("Output path for peak genomve coverage plot:", covplotPath, "\n")
##
cat("Generating genome-wide peak coverage plot...", "\n")
svg(file = covplotPath) # set the filepath for saving the svg figure
##
weightname <- names(narrowPeaks@elementMetadata@listData[2])
covplot(narrowPeaks, weightCol=weightname)
## Turn off svg device 
dev.off()


## Profile of peaks in TSS regions
TSSprofilePath <- gsub("peakgenomecov", "peakTSSprofile", covplotPath)
cat("Output path for peak TSS profile plot:", TSSprofilePath, "\n")
## One step function to generate TSS heatmap from a BED file
cat("Generating peak TSS profile plot...", "\n")
svg(file = TSSprofilePath) # set the filepath for saving the svg figure
##
peakHeatmap(bedFile, TxDb=txdb, upstream=3000, downstream=3000, color="red")
## Turn off svg device 
dev.off()


## Average profile of ChIP peaks binding to TSS region
AvgPeakProfileTSS <- gsub("peakgenomecov", "avgPeakTSSprofile", covplotPath)
cat("Output path for average peak TSS profile plot:", AvgPeakProfileTSS, "\n")
## One step function of above from a BED file
cat("Generating average peak TSS profile plot...", "\n")
svg(file = AvgPeakProfileTSS) # set the filepath for saving the svg figure
##
plotAvgProf2(bedFile, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
## Turn off svg device 
dev.off()


## With confidence interval estimate and bootstrap methods
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)


## Peak annotations
cat("Generating peak annotations...", "\n")
peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
## Plot it
annoPlot1 <- gsub("peakgenomecov", "annoplot1", covplotPath)
annoPlot2 <- gsub("peakgenomecov", "annoplot2", covplotPath)
annoPlot3 <- gsub("peakgenomecov", "annoplot3", covplotPath)
annoPlot4 <- gsub("peakgenomecov", "annoplot4", covplotPath)
#annoPlot5 <- gsub("peakgenomecov", "annoplot1", covplotPath)
##
svg(file = annoPlot1) # set the filepath for saving the svg figure
plotAnnoPie(peakAnno)
dev.off()
##
svg(file = annoPlot2) # set the filepath for saving the svg figure
plotAnnoBar(peakAnno)
dev.off()
##
svg(file = annoPlot3) # set the filepath for saving the svg figure
vennpie(peakAnno)
dev.off()
##
svg(file = annoPlot4) # set the filepath for saving the svg figure
upsetplot(peakAnno)
dev.off()
##
#svg(file = annoPlot5) # set the filepath for saving the svg figure
#upsetplot(peakAnno, vennpie=TRUE)


## Finish up
file.create(outputPath)


