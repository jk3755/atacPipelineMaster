#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  fastaPath <- snakemake@input[[1]]
  refBEDpath <- snakemake@input[[2]]
  bamPath <- snakemake@input[[3]]
  functionSourcePath <- snakemake@input[[4]]
  #
  modelOutPath <- snakemake@output[[1]]
  biasedPlotPath <- snakemake@output[[2]]
  correctedPlotPath <- snakemake@output[[3]]
  #
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Report ####
  cat("Generating seqbias model:", "\n")
  cat("Input fasta:", fastaPath, "\n")
  cat("Input BED interval file:", refBEDpath, "\n")
  cat("Input bam:", bamPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output path for seqbias model:", modelOutPath, "\n")
  cat("Output path for biased frequency plot:", biasedPlotPath, "\n")
  cat("Output path for corrected frequency plot:", correctedPlotPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", modelOutPath, "\n")
  if (file.exists(modelOutPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    cat("Output file not found. Processing", "\n")
  
  #### Load libraries ####
  cat("Loading libraries", "\n")
  suppressMessages(library(seqbias))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(ggplot2))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(genomation))
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  ####
  cat("Generating the seqbias model", "\n")
  seqbiasModel <- generateTn5BiasCorrectionModel(refFastaPath = fastaPath, bamPath = bamPath, refBEDpath = refBEDpath, biasedPlotPath = biasedPlotPath, correctedPlotPath = correctedPlotPath, biasModelPath = modelOutPath)
  
  ####
  cat("Finished", "\n")
  }
}, finally = {
})
