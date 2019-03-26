## See below link for reference information related to ChIPseeker package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Install libraries, if necessary
#install.packages("BiocManager")
#library(BiocManager)
#iocManager::install("ChIPseeker")
#BiocManager::install("genomation")
#BiocManager::install("GenomicRanges")
#BiocManager::install("clusterProfiler")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")

##
library(ChIPseeker)
library(genomation)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
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



