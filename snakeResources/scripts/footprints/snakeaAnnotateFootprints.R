##

## Load libraries #####################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
##
BiocManager::install("ChIPseeker")
##
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)


##
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##
x <- annotatePeak(
                  peak = test,
                  tssRegion = c(-3000, 3000),
                  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                  level = "transcript",
                  assignGenomicAnnotation = TRUE,
                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                  annoDb = NULL,
                  addFlankGeneInfo = FALSE,
                  flankDistance = 5000,
                  sameStrand = FALSE,
                  ignoreOverlap = FALSE,
                  ignoreUpstream = FALSE,ignoreDownstream = FALSE,
                  overlap = "TSS",
                  verbose = TRUE)







