## Set snakemake variables
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
dirPath <- snakemake@wildcards[["path"]]

## Report
cat("Bam filepath:", bamPath, "\n")
cat("Bai filepath:", baiPath, "\n")
cat("Output filepath:", outPath, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")
cat("Directory path:", dirPath, "\n")

## Load libraries
cat("Loading libraries", "\n")
suppressMessages(library(Rsamtools))

## Count reads in bam file
cat("Counting total sample reads", "\n")
bamFile <- BamFile(bamPath)
sampleTotalReads <- countBam(bamFile)
cat("Found", sampleTotalReads[[6]], "total reads in sample library", "\n")

## Save the data
cat("Saving data", "\n")
save(sampleTotalReads, file = outPath)