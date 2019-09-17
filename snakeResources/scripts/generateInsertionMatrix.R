#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  bamPath <- snakemake@input[[1]]
  baiPath <- snakemake@input[[2]]
  peakPath <- snakemake@input[[3]]
  insertionMatrixFilepath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  functionSourcePath <- snakemake@input[[4]]
  dirPath <- snakemake@wildcards[["path"]]
  
  #### Report ####
  cat("Generating insertion matrix with the following parameters:", "\n")
  cat("Bam file:", bamPath, "\n")
  cat("Bai file:", baiPath, "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Directory path:", dirPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Insertion matrix output path:", insertionMatrixFilepath, "\n")
  
  #### Perform a filecheck ####
  
  cat("Output filepath for insertion matrix:", insertionMatrixFilepath, "\n")
  
  if (file.exists(insertionMatrixFilepath) == TRUE){
    
    cat("Insertion matrix data file already exists, skipping operation", "\n")
    
  } else {
    
    cat("Insertion matrix data file not found, generating now", "\n")
    
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
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### Retrieve the binding sites for the current gene ####
    cat("Retrieving binding sites", "\n")
    tempAllSites <- getAllBindingSites(gene = geneName)
    
    #### Subset the binding sites with the peaks
    cat("Subsetting sites to peak regions", "\n")
    peaksGR <- importBED(peakPath)
    allSites <- subsetByOverlaps(tempAllSites, peaksGR)
    numSites <- length(allSites@ranges)
    cat("Identified", numSites, "total binding sites", "\n")
    
    if (numSites == 0){
      
    } else {
      
      #### Set scope of analysis ####
      maxWidth <- max(allSites@ranges@width)
      scope <- paste0("chr", c(1:22, "X", "Y"))
      cat("scope of analysis:", scope, "\n")
      
      #### Determine scope for current binding sites ####
      currentScope <- scope[which(scope %in% allSites@seqnames@values)]
      cat("scope of current binding sites:", currentScope, "\n")
      
      #### Generate insertion matrix for all sites ####
      for (item in currentScope)
      {
        com <- paste0("tempSites <- allSites[seqnames(allSites) == '", item, "']")
        eval(parse(text = com))
        com <- paste0(item, "insMatrix <- generateInsertionMatrixByChr(bamPath, tempSites, maxWidth,'", item, "')")
        eval(parse(text = com))
      }
      
      #### Merge the individual insertion matrices and save ####
      currentMatrixNames <- as.character(paste0(currentScope, "insMatrix"))
      currentMatrixNamesString <- paste(currentMatrixNames, collapse = ",")
      com <- paste0("insertionMatrix <- rbind(", currentMatrixNamesString, ")")
      eval(parse(text = com))
      insertionMatrixData <- list()
      insertionMatrixData$bindingSites <- allSites
      insertionMatrixData$insertionMatrix <- insertionMatrix

      #### Save the insertion matrix and input sites ####
      cat("Saving insertion matrix file", "\n")
      save(insertionMatrixData, file = insertionMatrixFilepath)
      
    }
  } # end filecheck

}, finally = {
  
})