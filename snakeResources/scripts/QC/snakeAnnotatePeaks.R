## See below link for reference information related to ChIPseeker package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

##
cat("Loading libraries", "\n")
if(!require(ChIPseeker)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ChIPseeker")}
if(!require(genomation)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("genomation")}
if(!require(GenomicRanges)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")}
if(!require(clusterProfiler)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("clusterProfiler")}
if(!require(org.Hs.eg.db)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")}
if(!require(ReactomePA)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ReactomePA")}
if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")}

## 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## Note that the .narrowPeak files are BED formatted
cat("Setting snakemake variables", "\n")
bedFile <- snakemake@input[[1]]
outputPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

## Import the peaks information from the .bed file
cat("Input peak file path:", bedFile, "\n")
narrowPeaks <- ChIPseeker::readPeakFile(bedFile)
## Subset to standard chromosomes
narrowPeaks <- keepStandardChromosomes(narrowPeaks, pruning.mode="coarse")


#### Coverage plots make plot of genome wide peak coverage
covplotPath <- gsub("annotations.done", "peak.genomecov.svg", outputPath)
covplotPath <- gsub("operations/preprocessing", "metrics", covplotPath)
cat("Output path for peak genomve coverage plot:", covplotPath, "\n")

##
cat("Generating genome-wide peak coverage plot", "\n")
svg(file = covplotPath) # set the filepath for saving the svg figure

##
weightname <- names(narrowPeaks@elementMetadata@listData[2])
covplot(narrowPeaks, weightCol = weightname)

## Turn off svg device 
dev.off()


#### Profile of peaks in TSS regions
TSSprofilePath <- gsub("peak.genomecov", "peak.TSSprofile", covplotPath)
cat("Output path for peak TSS profile plot:", TSSprofilePath, "\n")

## One step function to generate TSS heatmap from a BED file
cat("Generating peak TSS profile plot", "\n")
svg(file = TSSprofilePath) # set the filepath for saving the svg figure

##
peakHeatmap(bedFile, TxDb=txdb, upstream=3000, downstream=3000, color="red")

## Turn off svg device 
dev.off()


#### Average profile of ChIP peaks binding to TSS region
AvgPeakProfileTSS <- gsub("peak.genomecov", "avgPeak.TSSprofile", covplotPath)
cat("Output path for average peak TSS profile plot:", AvgPeakProfileTSS, "\n")

## One step function of above from a BED file
cat("Generating average peak TSS profile plot", "\n")

##
svg(file = AvgPeakProfileTSS) # set the filepath for saving the svg figure

##
plotAvgProf2(bedFile, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

## Turn off svg device 
dev.off()


#### Peak annotations
cat("Generating peak annotations", "\n")

##
peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

## Plot it
annoPlot1 <- gsub("peak.genomecov", "annoplot1", covplotPath)
annoPlot2 <- gsub("peak.genomecov", "annoplot2", covplotPath)
annoPlot3 <- gsub("peak.genomecov", "annoplot3", covplotPath)
annoPlot4 <- gsub("peak.genomecov", "annoplot4", covplotPath)

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


#### Finish up
file.create(outputPath)
