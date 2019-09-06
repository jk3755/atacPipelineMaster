## Set snakemake variables
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outputSVG <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

## Report
cat("Bam filepath:", bamPath, "\n")
cat("Bai filepath:", baiPath, "\n")
cat("Output filepath:", outputSVG, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")

## Load libraries
cat("Loading libraries", "\n")
suppressMessages(library(ATACseqQC))

##
cat("Generating fragment size distribution", "\n")
sampleLabel <- paste0(sampleName, "-", sampleRep)
svg(file = outputSVG)
fragSizeDist(
  bamFiles = bamPath,
  bamFiles.labels = sampleLabel,
  index = baiPath)
dev.off()