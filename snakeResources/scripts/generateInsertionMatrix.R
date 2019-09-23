#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  bamPath <- snakemake@input[[1]]
  baiPath <- snakemake@input[[2]]
  functionSourcePath <- snakemake@input[[3]]
  bindingSitesPath <- snakemake@input[[4]]
  snakeTouchPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Report ####
  cat("Generating insertion matrix with the following parameters:", "\n")
  cat("Bam file:", bamPath, "\n")
  cat("Bai file:", baiPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Filepath for loading binding sites:", bindingSitesPath, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", snakeTouchPath, "\n")
  if (file.exists(snakeTouchPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    cat("Output file not found. Processing", "\n")
    
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
    suppressMessages(library(stringr))
    genome <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    #### Generate the output directory path ####
    outputDirectory <- str_replace(snakeTouchPath, "(?<=insertion_matrix/).*", paste0(geneName, "/"))
    cat("Generating output directory at path:", outputDirectory, "\n")
    dir.create(path = outputDirectory, showWarnings = FALSE)
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### Retrieve the binding sites for the current gene ####
    cat("Loading binding sites", "\n")
    load(bindingSitesPath)
    
    numSites <- length(bindingSites@ranges)
    cat("Identified", numSites, "total binding sites", "\n")
    
    if (numSites == 0){
      cat("No binding sites identified, exiting", "\n")
    } else {
      
      #### Determine scope for current binding sites ####
      maxWidth <- max(bindingSites@ranges@width)
      scope <- paste0("chr", c(1:22, "X", "Y"))
      cat("scope of analysis:", scope, "\n")
      currentScope <- scope[which(scope %in% bindingSites@seqnames@values)]
      cat("scope of current binding sites:", currentScope, "\n")
      
      #### Generate insertion matrix for all sites ####
      for (item in currentScope)
      {
        ##
        cat("Processing insertion matrix for", item, "\n")
        
        ## Subset the binding sites
        cat("Subsetting binding sites", "\n")
        com <- paste0("currentSites <- bindingSites[seqnames(bindingSites) == '", item, "']")
        eval(parse(text = com))
        
        ## Generate the matrix
        cat("Generating insertion matrix", "\n")
        com <- paste0("insMatrix <- generateInsertionMatrixByChr(bamPath, currentSites, maxWidth,'", item, "')")
        eval(parse(text = com))
        
        ## Make the list
        cat("Storing data", "\n")
        insertionMatrixData <- list()
        insertionMatrixData$bindingSites <- currentSites
        insertionMatrixData$insertionMatrix <- insMatrix
        
        ## Save the file
        matrixOutPath <- paste0(outputDirectory, item, ".", geneName, ".insertionMatrix.RData")
        cat("Saving insertion matrix at path:", matrixOutPath, "\n")
        save(insertionMatrixData, file = matrixOutPath)
      }
      
      #### Touch the snakemake file ####
      cat("Finished, touching snakemake file", "\n")
      file.create(snakeTouchPath, showWarnings = FALSE)
      
    }
  }
  #### Touch the snakemake file ####
  cat("Finished, touching snakemake file", "\n")
  file.create(snakeTouchPath, showWarnings = FALSE)
  
}, finally = {
  
  #### Touch the snakemake file ####
  cat("Finished, touching snakemake file", "\n")
  file.create(snakeTouchPath, showWarnings = FALSE)
  
})