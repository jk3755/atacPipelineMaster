##
library(VariantAnnotation)

##
rrbsPath <- "C:\\Users\\jsk33\\Desktop\\rrbs\\SRR8633497_1.fastq_hg38.fa.methylcytosines.vcf"

##
rrbs <- readVcf(file = rrbsPath, genome = "hg38")

## Trim the rrbs data to standard chromosomes only
scope <- paste0("chr", c(1:22, "X", "Y"))
rrbs <- keepStandardChromosomes(rrbs, pruning.mode="coarse")
rrbs <- keepSeqlevels(rrbs, scope, pruning.mode="coarse")
rrbs <- trim(rrbs, use.names = TRUE)
numSites <- length(rrbs)


## VCF methylation calls by bicycle
## CHG - methylation call for a cytosine in a CG context
## CHH - methylation call for a cytosine in a non-CG context


load("C:/Users/jsk33/Desktop/H508A-WT-02.ODC1.rawFootprintData.Rdata")