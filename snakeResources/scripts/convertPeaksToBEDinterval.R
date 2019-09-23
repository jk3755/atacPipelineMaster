#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  peaksPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  BEDoutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Report ####
  cat("Generating seqbias model:", "\n")
  cat("Input peaks:", peaksPath, "\n")
  cat("Output path for BED:", BEDoutPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  
  #### Load libraries ####
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(genomation))
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  #### Load BED intervals ####
  cat("Loading peaks file", "\n")
  peaksGR <- importBED(peaksPath)
  
  #### Load BED intervals ####
  cat("Writing peaks to BED interval file", "\n")
  writeGRangesToBEDwithStrand(peaksGR, BEDoutPath)
  
  ####
  cat("Finished", "\n")
  
}, finally = {
})
