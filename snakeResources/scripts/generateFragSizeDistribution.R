
##
cat("Loading libraries", "\n")
suppressMessages(library(ATACseqQC))

##
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outputSVG <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

##
sampleLabel <- paste0(sampleName, "-", sampleRep)

##
cat("Output filepath for fragment size distribution: ", outputSVG, "\n")
svg(file = outputSVG) # set the filepath for saving the svg figure

## Calculate the fragment size distribution and save to svg
fragSizeDist(
  bamFiles = bamPath,
  bamFiles.labels = sampleLabel,
  index = baiPath
)

## Turn off svg device 
dev.off()

