#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Memory statistics ####
  cat("Setting gc options", "\n")
  gcinfo(TRUE)
  gc()
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  bamPath <- snakemake@input[[1]]
  baiPath <- snakemake@input[[2]]
  peakPath <- snakemake@input[[3]]
  snakemakeTouchPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  functionSourcePath <- snakemake@input[[4]]
  dirPath <- snakemake@wildcards[["path"]]
  
  #### Report ####
  cat("Spooling footprint analysis", "\n")
  cat("Bam file:", bamPath, "\n")
  cat("Bai file:", baiPath, "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Directory path:", dirPath, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Snakemake touchfile path:", snakemakeTouchPath, "\n")
  
  #### Perform a filecheck ####
  footprintDataFilepath <- paste0(dirPath, "footprints/samplemergedpeaks/raw/", sampleName, ".rep", sampleRep, ".ref", refGenome, ".", geneName, ".sampleMergedPeaks.FPdata.Rdata")
  cat("Output filepath for footprint data:", footprintDataFilepath, "\n")
  
  if (file.exists(footprintDataFilepath) == TRUE){
    cat("Footprint data file already exists, skipping operation", "\n")
  } else {
    cat("Footprint data file not found, generating", "\n")
    
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
      
      maxWidth <- max(allSites@ranges@width)
      
      #### Set scope of analysis ####
      scope <- paste0("chr", c(1:22, "X", "Y"))
      cat("scope of analysis:", scope, "\n")
      
      #### Determine scope for current binding sites ####
      currentScope <- scope[which(scope %in% allSites@seqnames@values)]
      cat("scope of current binding sites:", currentScope, "\n")
      
      #### Generate insertion matrix for each chromosome one at a time and convert to basic footprint statistics ####
      ## Doing it this way allows the script to keep the memory usage low
      for (item in currentScope){
        cat("Generating insertion matrix for", item, "\n")
        com <- paste0("tempSites <- allSites[seqnames(allSites) == '", item, "']")
        eval(parse(text = com))
        com <- paste0("insMatrix <- generateInsertionMatrixByChr(bamPath,tempSites,maxWidth,'", item, "')")
        eval(parse(text = com))
        com <- paste0(item, "SiteStatistics <- calculateBasicFootprintStatistics(insMatrix,tempSites,'", item, "')")
        eval(parse(text = com))
        rm(insMatrix)
      }
      
      #### Clear memory ####
      rm(allSites)
      gc()
      
      #### Merge all of the site statistics tables ####
      cat("Merging site statistics tables", "\n")
      currentTables <- as.character(paste0(currentScope, "SiteStatistics"))
      currentTablesStr <- paste(currentTables, collapse = ",")
      com <- paste0("footprintSiteStatistics <- rbind(", currentTablesStr, ")")
      eval(parse(text = com))
      
      cat("Saving footprint data file", "\n")
      save(footprintSiteStatistics, file = footprintDataFilepath)
    }
  } # end filecheck
  
  cat("Finished, touching snakemake flag file", "\n")
  file.create(snakemakeTouchPath)
  
}, finally = {
  
  cat("Finished, touching snakemake flag file", "\n")
  file.create(snakemakeTouchPath)
})