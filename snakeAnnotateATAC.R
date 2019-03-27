## See below link for reference information related to ChIPseeker package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Install libraries, if necessary
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("ChIPseeker")
#BiocManager::install("genomation")
#BiocManager::install("GenomicRanges")
#BiocManager::install("clusterProfiler")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")

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
svg(file = outputSVG) # set the filepath for saving the svg figure
##
covplot(narrowPeaks, weightCol="X486")
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
annoPlot2 <- gsub("peakgenomecov", "annoplot1", covplotPath)
annoPlot3 <- gsub("peakgenomecov", "annoplot1", covplotPath)
annoPlot4 <- gsub("peakgenomecov", "annoplot1", covplotPath)
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

## Distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")

## Functional enrichment analysis
pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)
#
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
#
dotplot(pathway2)


## Profile of several peak data sets binding to TSS
## promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
## tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
##
## to speed up the compilation of this vigenette, we load a precaculated tagMatrixList
data("tagMatrixList")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

##
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")


##
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)




## Peak annotation comparisons
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

## Functional profiles comparison
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")



## Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)


## Statistical testing of peak overlaps
# shuffle genome coordination
p <- GRanges(seqnames=c("chr1", "chr3"),
             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)
## Peak overlap enrichment analysis
enrichPeakOverlap(queryPeak     = files[[5]],
                  targetPeak    = unlist(files[1:4]),
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)


## Data mining with ChIP seq data deposited in GEO
getGEOspecies()
getGEOgenomeVersion()
##
hg38 <- getGEOInfo(genome="hg38", simplify=TRUE)
head(hg38)

## Download GEO ChIP datasets
downloadGEObedFiles(genome="hg19", destDir="hg19")
#
gsm <- hg19$gsm[sample(nrow(hg19), 10)]
downloadGSMbedFiles(gsm, destDir="hg19")



