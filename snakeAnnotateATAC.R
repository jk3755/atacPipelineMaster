

##
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("ChIPseeker")
BiocManager::install("genomation")
BiocManager::install("GenomicRanges")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")


##
library(ChIPseeker)
library(genomation)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


## Data import
cat("Importing peaks file...", "\n")
bedFile <- "C:\\Users\\jsk33\\Documents\\atac\\atac1\\mdst8\\wt01\\peaks\\macs2\\merged\\MDST8-WT-01-merged_local_normalization_peaks.narrowPeak"
narrowPeaks <- ChIPseeker::readPeakFile(bedFile)
## Subset to standard xsomes
narrowPeaks <- keepStandardChromosomes(narrowPeaks, pruning.mode="coarse")


## Coverage plots make plot of genome wide peak coverage
covplot(narrowPeaks, weightCol="X9.11849")


## More detailed cov plots
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))


## Profile of peaks in TSS regions
#promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#tagMatrix <- getTagMatrix(bedFile, windows=promoter)
## heatmap of binding to TSS regions
#tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

## One step function to generate TSS heatmap from a BED file
peakHeatmap(bedFile, TxDb=txdb, upstream=3000, downstream=3000, color="red")


## Average profile of ChIP peaks binding to TSS region
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


## One step function of above from a BED file
plotAvgProf2(bedFile, TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

## With confidence interval estimate and bootstrap methods
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)


## Peak annotation
peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
## Plot it
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)






