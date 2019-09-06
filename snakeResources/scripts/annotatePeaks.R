## Set snakemake variables
cat("Setting snakemake variables", "\n")
bedFile <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

## Report
cat("BED filepath:", bedFile, "\n")
cat("Output filepath:", outPath, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")

## Load libraries
cat("Loading libraries", "\n")
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicRanges))
suppressMessages(library(org.Hs.eg.db))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Import the peaks information from the .bed file
cat("Importing peaks from BED file", "\n")
narrowPeaks <- ChIPseeker::readPeakFile(bedFile)
narrowPeaks <- keepStandardChromosomes(narrowPeaks)
narrowPeaks <- trim(narrowPeaks)

## Coverage plots make plot of genome wide peak coverage
cat("Generating genome coverage plots", "\n")
covplotPath <- gsub("annotations.done", "peak.genomecov.svg", outPath)
covplotPath <- gsub("operations/metrics", "metrics", covplotPath)
cat("Output path for peak genomve coverage plot:", covplotPath, "\n")
svg(file = covplotPath)
weightname <- names(narrowPeaks@elementMetadata@listData[2])
covplot(narrowPeaks, weightCol = weightname)
dev.off()

## Profile of peaks in TSS regions
cat("Generating TSS peak profile", "\n")
TSSprofilePath <- gsub("peak.genomecov", "peak.TSSprofile", covplotPath)
cat("Output path for peak TSS profile plot:", TSSprofilePath, "\n")
svg(file = TSSprofilePath)
peakHeatmap(bedFile, TxDb=txdb, upstream=3000, downstream=3000, color="red")
dev.off()

## Average profile of ChIP peaks binding to TSS region
cat("Generating average TSS peak profile", "\n")
AvgPeakProfileTSS <- gsub("peak.genomecov", "avgPeak.TSSprofile", covplotPath)
cat("Output path for average peak TSS profile plot:", AvgPeakProfileTSS, "\n")
cat("Generating average peak TSS profile plot", "\n")
svg(file = AvgPeakProfileTSS)
plotAvgProf2(bedFile, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

## Peak annotations
cat("Generating peak annotations", "\n")
peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
annoPlot1 <- gsub("peak.genomecov", "annoplot1", covplotPath)
annoPlot2 <- gsub("peak.genomecov", "annoplot2", covplotPath)
annoPlot3 <- gsub("peak.genomecov", "annoplot3", covplotPath)
annoPlot4 <- gsub("peak.genomecov", "annoplot4", covplotPath)
svg(file = annoPlot1)
plotAnnoPie(peakAnno)
dev.off()
svg(file = annoPlot2)
plotAnnoBar(peakAnno)
dev.off()
svg(file = annoPlot3)
vennpie(peakAnno)
dev.off()
svg(file = annoPlot4)
upsetplot(peakAnno)
dev.off()

## Finish up
cat("Touching snakemake flag file", "\n")
file.create(outPath)