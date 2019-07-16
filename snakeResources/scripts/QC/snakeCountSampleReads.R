## Set snakemake variables
cat("Setting snakemake variables...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
dirPath <- snakemake@wildcards[["path"]]

##
# cat("Counting total reads in library for sample:", sampleName, "\n")
# if(!require(Rsamtools)){
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("Rsamtools")}

##
bamFile <- BamFile(bamPath)
sampleTotalReads <- countBam(bamFile)
sampleTotalReads <- sampleTotalReads[[6]]
cat("Found", sampleTotalReads, "total reads in sample library", "\n")
      
## Save the data
save(sampleTotalReads, file = outPath)