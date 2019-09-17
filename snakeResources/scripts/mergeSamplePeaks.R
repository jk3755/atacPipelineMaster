
#### Determine the number of input files ####
lastInput <- length(snakemake@input)
numInputs <- (lastInput - 1)
cat("Processing", numInputs, "input files", "\n")

#### Set snakemake variables ####
functionSourcePath <- snakemake@input[[lastInput]]
outPath <- snakemake@output[[1]]
for (a in 1:numInputs){
  com <- paste0("inputPath", a, " <- snakemake@input[[", a, "]]")
  eval(parse(text = com))}

#### Load libraries ####
cat("Loading libraries", "\n")
suppressMessages(library(genomation))
suppressMessages(library(GenomicRanges))

#### Source functions ####
cat("Loading functions from path:", functionSourcePath, "\n")
source(functionSourcePath)

#### Load peaks from .bed file to GRanges object ####
cat("Loading peaks from .bed files as GRanges", "\n")
for (b in 1:numInputs){
  com <- paste0("samplePeaks", b, " <- importBED(inputPath", b, ")")
  eval(parse(text = com))}

#### Transfer GRanges to a list to use the union function ####
cat("Storing all GRanges object in a list", "\n")
peaksList <- list()
for (c in 1:numInputs){
  com <- paste0("peaksList[[", c, "]] <- samplePeaks", c)
  eval(parse(text = com))}

#### Merge the genomic ranges by using union function ####
cat("Finding the union of all GRanges objects", "\n")
peaksUnion <- getGRangesUnion(peaksList)

#### Write the output file ####
cat("Writing merged peaks data to file at:", outPath, "\n")
writeGRangesToBED(peaksUnion, outPath)
