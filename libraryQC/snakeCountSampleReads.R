## Set snakemake variables
cat("Setting snakemake variables...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
dirPath <- snakemake@wildcards[["path"]]

##
cat("Counting total reads in library for sample:", sampleName, "\n")
suppressMessages(library(Rsamtools))

##
bamFile <- BamFile(bamPath)
sampleTotalReads <- countBam(bamFile)
cat("Found", sampleTotalReads, "total reads in sample library", "\n")
      
## Save the data
cat("Saving finished raw footprint data for", geneName, "\n")
save(sampleTotalReads, file = outPath)