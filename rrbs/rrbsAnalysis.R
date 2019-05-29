##
library(VariantAnnotation)
library(GenomicRanges)

##
rrbsPath <- "C:\\Users\\jsk33\\Desktop\\SRR8633497_1.fastq_hg38.fa.methylcytosines.vcf"

## VCF methylation calls by bicycle
## CHG - methylation call for a cytosine in a CG context
## CHH - methylation call for a cytosine in a non-CG context
## CG - non-methylated?
rrbs <- readVcf(file = rrbsPath, genome = "hg38")

## Trim the rrbs data to standard chromosomes only
scope <- paste0("chr", c(1:22, "X", "Y"))
rrbs <- keepStandardChromosomes(rrbs, pruning.mode="coarse")
rrbs <- keepSeqlevels(rrbs, scope, pruning.mode="coarse")
rrbs <- trim(rrbs, use.names = TRUE)
numSites <- length(rrbs)

##
load("C:/Users/jsk33/Desktop/H508A-WT-02.ODC1.rawFootprintData.Rdata")

## FPs
grFP <- footprintData[["motif1"]][["genomeSites"]]


## Overlaps
methOverlaps <- findOverlaps(grFP, rrbs)
methIdx <- methOverlaps@to
meth2 <- rrbs[methIdx]

##
x <- which(meth2@rowRanges@ranges@NAMES == "CHG")
y <- which(meth2@rowRanges@ranges@NAMES == "CHH")
z <- which(meth2@rowRanges@ranges@NAMES == "CG")
