#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  inputString <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  dataOutPath <- snakemake@output[[1]]
  sampleDirectory <- snakemake@wildcards[["path"]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Report ####
  cat("Spooling raw footprint analysis", "\n")
  cat("Input string:", inputString, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for aggregated footprint data:", dataOutPath, "\n")
  cat("Sample directory:", sampleDirectory, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Reference genome used:", refGenome, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", dataOutPath, "\n")
  if (file.exists(dataOutPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    
    #### Detect processing type ####
    suppressMessages(library(stringr))
    sampleSpecific <- FALSE
    sampleMerged <- FALSE
    sampleSpecific <- str_detect(inputString, "sample_specific")
    sampleMerged <- str_detect(inputString, "sample_merged")
    ##
    if (sampleSpecific == TRUE){
      cat("Processing type detected as sample specific", "\n")
      inputDirectory <- paste0(sampleDirectory, "footprints/sample_specific/raw/")
    } else if (sampleMerged == TRUE){
      cat("Processing type detected as sample merged", "\n")
      inputDirectory <- paste0(sampleDirectory, "footprints/sample_merged/raw/")
    }
    
    ####
    cat("Input directory at path:", inputDirectory, "\n")
    inputFileList <- list.files(path = inputDirectory, full.names = TRUE)
    numInputFiles <- length(inputFileList)
    cat("Identified", numInputFiles, "input files", "\n")
    
    #### Load libraries ####
    cat("Loading libraries", "\n")
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(Biostrings))
    suppressMessages(library(MotifDb))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(Rsamtools))
    suppressMessages(library(GenomicAlignments))
    suppressMessages(library(ChIPseeker))
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    suppressMessages(library(genomation))
    genome <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    ####
    cat(inputFileList, "\n")
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### List for aggregating data ####
    aggregatedFootprints <- list()
    
    #### Analyze the footprints ####
    for (a in 1:numInputFiles){
      
      ## Get the gene name
      inputFile <- inputFileList[a]
      geneName <- str_extract(inputFile, "(?<=refhg38.)[\\w|\\d]*")
      cat("Extracted gene name:", geneName, "\n")
      
      ## Load footprint data
      cat("Loading data from path:", inputFile, "\n")
      load(file = inputFile)
      
      ## Transfer the data to the aggregate list
      com <- paste0("aggregatedFootprints$", geneName, " <- footprintSiteStatistics")
    }
    
    #### Save the aggregated footprint data ####
    cat("Saving aggregated footprint data", "\n")
    save(aggregatedFootprints, file = dataOutPath)
  }
}, finally = {
})
