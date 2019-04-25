
## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)

## Disable scientific notation in variables and turn off warnings globally
options(scipen = 999)
options(warn = -1)

## Set snakemake variables
cat("Setting snakemake variables", "\n")
inputPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]
bamCopy <- snakemake@wildcards[["bamcopy"]]

## Set the filepaths and perform a filecheck
cat("Setting filepaths and checking for output file", "\n")
inputPath <- gsub("operations/parse", "data/parsed", inputPath)
inputPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", inputPath, perl = TRUE)
dataOutPath <- gsub("parsed", "processed", inputPath)

tryCatch({

##
if (file.exists(dataOutPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  ## Load libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(rlist))
  
  ## Load the current parsed footprintData object
  cat("Loading parsed footprint file", "\n")
  load(inputPath)
  ## To avoid errors, clear the list of any empty sub-lists first
  cat("Cleaning list", "\n")
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  ## Also remove any sublists that do not have the parsed data object stored
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 15L, TRUE)
  numMotifs <- length(footprintData)
  cat("Found", numMotifs, "motifs", "\n")
  
  ## Error handling
  if (numMotifs == 0){
    
    cat("No motifs found. Exiting", "\n")
    
  } else if (numMotifs == 1){

  ## Because some motifs may be removed due to errors, pull the motif names to be used in downstream commands (can't just go sequentially)
  motifNames <- names(footprintData)
  
  #### MERGE AND DEDUPLICATE ALL BINDING SITES ####
  ## At this point, the footprints have already been parsed into bound and unbound sites
  ## Transfer data for both to the new storage object, for downstream analysis
  # If there is only one motif available, there is no need to merge and deduplicate the identified sites
    cat("Processing 1 motif", "\n")
    cat("Setting peak sites", "\n")
    com <- paste0("peakSites <- footprintData[['", motifNames[1], "']][['peakSites']]")
    eval(parse(text = com))
    cat("Setting peak insertion matrix", "\n")
    com <- paste0("peakInsertionMatrix <- footprintData[['", motifNames[1], "']][['insMatrix']]")
    eval(parse(text = com))
    cat("Setting footprint metrics", "\n")
    com <- paste0("rawPeakFootprintMetrics <- footprintData[['", motifNames[1], "']][['rawFootprintMetrics']]")
    eval(parse(text = com))
    numPeakSites <- length(peakSites)
    cat("Number of peak sites", numPeakSites, "\n")
    ## Split the sites into bound and unbound as determined by null model with bonferroni correction
    cat("Finding bound site overlaps", "\n")
    com <- paste0("boundSiteOverlaps <- findOverlaps(footprintData[['", motifNames[1], "']][['parseData']][['bfSites']], footprintData[['", motifNames[1], "']][['peakSites']])")
    eval(parse(text = com))
    cat("Setting index for bound sites", "\n")
    boundSiteIndex <- boundSiteOverlaps@to
    cat("Subsetting bound sites", "\n")
    boundSites <- peakSites[boundSiteIndex]
    boundSitesInsertionMatrix <- peakInsertionMatrix[boundSiteIndex,]
    boundSitesMetrics <- rawPeakFootprintMetrics[boundSiteIndex,]
    numBoundSites <- length(boundSites)
    cat("Subsetting unbound sites", "\n")
    unboundSites <- peakSites[-boundSiteIndex]
    unboundSitesInsertionMatrix <- peakInsertionMatrix[-boundSiteIndex,]
    unboundSitesMetrics <- rawPeakFootprintMetrics[-boundSiteIndex,]
    numUnboundSites <- length(unboundSites)
    
    #### Calculate footprint characteristics on merged data ####
    ## Calculate the 10% trimmed mean of all insertions in the motif sites
    peakMotifSignal <- mean(rawPeakFootprintMetrics[,3], trim = 0.10)
    boundMotifSignal <- mean(boundSitesMetrics[,3], trim = 0.10)
    unboundMotifSignal <- mean(unboundSitesMetrics[,3], trim = 0.10)
    
    ## Calculate the mean of all insertions in the flank region
    peakFlankSignal <- mean(rawPeakFootprintMetrics[,2], trim = 0.10)
    boundFlankSignal <- mean(boundSitesMetrics[,2], trim = 0.10)
    unboundFlankSignal <- mean(unboundSitesMetrics[,2], trim = 0.10)
    
    ## Calculate the mean of background insertions
    peakBackgroundSignal <- mean(rawPeakFootprintMetrics[,1], trim = 0.10)
    boundBackgroundSignal <- mean(boundSitesMetrics[,1], trim = 0.10)
    unboundBackgroundSignal <- mean(unboundSitesMetrics[,1], trim = 0.10)
    
    ## Calculate flanking accessibility (log2 fold change between flank and background)
    peak.log2Flank <- log2(peakFlankSignal / peakBackgroundSignal)
    bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
    unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
    
    ## Calculate footprint depth (log2 fold change between flank and background)
    peak.log2Depth <- log2(peakMotifSignal / peakFlankSignal)
    bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
    unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
    
    #### TRANSFER DATA TO STORAGE OBJECT ####
    ## Initialize a new list object to store the processed data
    processedFootprintData <- list()
    ##
    processedFootprintData$"geneName" <- footprintData[["motif1"]][["geneName"]]
    processedFootprintData$"numMotifs" <- numMotifs
    processedFootprintData$"numPeakSites" <- numPeakSites
    processedFootprintData$"numBoundSites" <- numBoundSites
    processedFootprintData$"numUnboundSites" <- numUnboundSites
    ##
    processedFootprintData$"peakSites" <- peakSites
    processedFootprintData$"rawPeakFootprintMetrics" <- rawPeakFootprintMetrics
    processedFootprintData$"peakMotifSignal" <- peakMotifSignal
    processedFootprintData$"peakFlankSignal" <- peakFlankSignal
    processedFootprintData$"peakBackgroundSignal" <- peakBackgroundSignal
    processedFootprintData$"peak.log2Flank" <- peak.log2Flank
    processedFootprintData$"peak.log2Depth" <- peak.log2Depth
    ##
    processedFootprintData$"boundSites" <- boundSites
    processedFootprintData$"boundSitesMetrics" <- boundSitesMetrics
    processedFootprintData$"boundMotifSignal" <- boundMotifSignal
    processedFootprintData$"boundFlankSignal" <- boundFlankSignal
    processedFootprintData$"boundBackgroundSignal" <- boundBackgroundSignal
    processedFootprintData$"bound.log2Flank" <- bound.log2Flank
    processedFootprintData$"bound.log2Depth" <- bound.log2Depth
    ##
    processedFootprintData$"unboundSites" <- unboundSites
    processedFootprintData$"unboundSitesMetrics" <- unboundSitesMetrics
    processedFootprintData$"unboundMotifSignal" <- unboundMotifSignal
    processedFootprintData$"unboundFlankSignal" <- unboundFlankSignal
    processedFootprintData$"unboundBackgroundSignal" <- unboundBackgroundSignal
    processedFootprintData$"unbound.log2Flank" <- unbound.log2Flank
    processedFootprintData$"unbound.log2Depth" <- unbound.log2Depth
    
    ## Save the data
    save(processedFootprintData, file = dataOutPath)
    
  } else {
    
    motifNames <- names(footprintData)
    ## The insertion matrices cannot be concatenated because they have different numbers of columns
    ## Initialize list objects here to store the insertion matrices separately (not concatenated)
    ## Pull the data for each individual motif and perform the bound/unbound split
    cat("Processing multiple motifs", "\n")
    ##
    for (z in 1:numMotifs){
      ## Pull the basic data
      cat("Pulling peak sites", "\n")
      com <- paste0("peakSites", z, " <- footprintData[['", motifNames[z], "']][['peakSites']]")
      eval(parse(text = com))
      cat("Pulling insertion matrix", "\n")
      com <- paste0("peakInsertionMatrix", z, " <- footprintData[['", motifNames[z], "']][['insMatrix']]")
      eval(parse(text = com))
      cat("Pulling footprint metrics", "\n")
      com <- paste0("rawPeakFootprintMetrics", z, " <- footprintData[['", motifNames[z], "']][['rawFootprintMetrics']]")
      eval(parse(text = com))
      ## Perform the bound/unbound split
      cat("Finding bound overlaps", "\n")
      com <- paste0("boundSiteOverlaps", z, " <- findOverlaps(footprintData[['", motifNames[z], "']][['parseData']][['bfSites']], footprintData[['", motifNames[z], "']][['peakSites']])")
      eval(parse(text = com))
      cat("Setting bound index", "\n")
      com <- paste0("boundSiteIndex", z, " <- boundSiteOverlaps", z, "@to")
      eval(parse(text = com))
      ## Bound sites
      cat("Setting bound site data", "\n")
      com <- paste0("boundSites", z, " <- peakSites", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("boundSitesInsertionMatrix", z, " <- peakInsertionMatrix", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("boundSitesMetrics", z, " <- rawPeakFootprintMetrics", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("numBoundSites", z, " <- length(boundSites", z, ")")
      eval(parse(text = com))
      ## Unbound sites
      cat("Setting unbound site data", "\n")
      com <- paste0("unboundSites", z, " <- peakSites", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("unboundSitesInsertionMatrix", z, " <- peakInsertionMatrix", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("unboundSitesMetrics", z, " <- rawPeakFootprintMetrics", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("numUnboundSites", z, " <- length(unboundSites", z, ")")
      eval(parse(text = com))
    } # end for (z in 1:numMotif)
    
    #### Perform the merging and deduplication ####
    for (b in 2:numMotifs){
      ## Find overlaps and generate selection indices
      com <- paste0("tempOverlapsPeaks <- findOverlaps(peakSites1, peakSites", b,")")
      eval(parse(text = com))
      com <- paste0("tempOverlapsBound <- findOverlaps(boundSites1, boundSites", b,")")
      eval(parse(text = com))
      com <- paste0("tempOverlapsUnbound <- findOverlaps(unboundSites1, unboundSites", b,")")
      eval(parse(text = com))
      
      #### MERGE THE PEAK SITES (all sites) ####
      ## If no overlaps are present, can just directly merge the two Granges
      ## Otherwise, if some overlaps are present, merge the Granges, but omit overlapping sites from second group
      if (length(tempOverlapsPeaks@from) == 0){
        com <- paste0("peakSites1 <- c(peakSites1, peakSites", b, ")")
        eval(parse(text = com))
        com <- paste0("rawPeakFootprintMetrics1 <- rbind(rawPeakFootprintMetrics1, rawPeakFootprintMetrics", b, ")")
        eval(parse(text = com))
      } else {
        mergeIdx <- tempOverlapsPeaks@to
        com <- paste0("peakSites1 <- c(peakSites1, peakSites", b, "[-mergeIdx])")
        eval(parse(text = com))
        com <- paste0("rawPeakFootprintMetrics1 <- rbind(rawPeakFootprintMetrics1, rawPeakFootprintMetrics", b, "[-mergeIdx,])")
        eval(parse(text = com))
      } # end if (length(overlaps@from) == 0)
      
      #### MERGE THE BOUND SITES ####
      if (length(tempOverlapsBound@from) == 0){
        com <- paste0("boundSites1 <- c(boundSites1, boundSites", b, ")")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics1 <- rbind(boundSitesMetrics1, boundSitesMetrics", b, ")")
        eval(parse(text = com))
      } else {
        mergeIdx <- tempOverlapsBound@to
        com <- paste0("boundSites1 <- c(boundSites1, boundSites", b, "[-mergeIdx])")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics1 <- rbind(boundSitesMetrics1, boundSitesMetrics", b, "[-mergeIdx,])")
        eval(parse(text = com))
      } # end if (length(tempOverlapsBound@from) == 0)
      
      #### MERGE THE UNBOUND SITES ####
      if (length(tempOverlapsUnbound@from) == 0){
        com <- paste0("unboundSites1 <- c(unboundSites1, unboundSites", b, ")")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics1 <- rbind(unboundSitesMetrics1, unboundSitesMetrics", b, ")")
        eval(parse(text = com))
      } else {
        mergeIdx <- tempOverlapsUnbound@to
        com <- paste0("unboundSites1 <- c(unboundSites1, unboundSites", b, "[-mergeIdx])")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics1 <- rbind(unboundSitesMetrics1, unboundSitesMetrics", b, "[-mergeIdx,])")
        eval(parse(text = com))
      } # end if (length(tempOverlapsUnbound@from) == 0)
    } # end for (b in 2:numMotifs)
    
    ## Because the code is written to use the first motif to merge everything into,
    ## transfer the data at this stage to make it consistent with numMotifs == 1
    peakSites <- peakSites1
    rawPeakFootprintMetrics <- rawPeakFootprintMetrics1
    numPeakSites <- length(peakSites)
    ##
    boundSites <- boundSites1
    boundSitesMetrics <- boundSitesMetrics1
    numBoundSites <- length(boundSites)
    ##
    unboundSites <- unboundSites1
    unboundSitesMetrics <- unboundSitesMetrics1
    numUnboundSites <- length(unboundSites)
    
  #### Calculate footprint characteristics on merged data ####
  ## Calculate the 10% trimmed mean of all insertions in the motif sites
  peakMotifSignal <- mean(rawPeakFootprintMetrics[,3], trim = 0.10)
  boundMotifSignal <- mean(boundSitesMetrics[,3], trim = 0.10)
  unboundMotifSignal <- mean(unboundSitesMetrics[,3], trim = 0.10)
  
  ## Calculate the mean of all insertions in the flank region
  peakFlankSignal <- mean(rawPeakFootprintMetrics[,2], trim = 0.10)
  boundFlankSignal <- mean(boundSitesMetrics[,2], trim = 0.10)
  unboundFlankSignal <- mean(unboundSitesMetrics[,2], trim = 0.10)
  
  ## Calculate the mean of background insertions
  peakBackgroundSignal <- mean(rawPeakFootprintMetrics[,1], trim = 0.10)
  boundBackgroundSignal <- mean(boundSitesMetrics[,1], trim = 0.10)
  unboundBackgroundSignal <- mean(unboundSitesMetrics[,1], trim = 0.10)
  
  ## Calculate flanking accessibility (log2 fold change between flank and background)
  peak.log2Flank <- log2(peakFlankSignal / peakBackgroundSignal)
  bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
  unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
  
  ## Calculate footprint depth (log2 fold change between flank and background)
  peak.log2Depth <- log2(peakMotifSignal / peakFlankSignal)
  bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
  unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
  
  #### TRANSFER DATA TO STORAGE OBJECT ####
  ## Initialize a new list object to store the processed data
  processedFootprintData <- list()
  ##
  processedFootprintData$"geneName" <- footprintData[["motif1"]][["geneName"]]
  processedFootprintData$"numMotifs" <- numMotifs
  processedFootprintData$"numPeakSites" <- numPeakSites
  processedFootprintData$"numBoundSites" <- numBoundSites
  processedFootprintData$"numUnboundSites" <- numUnboundSites
  ##
  processedFootprintData$"peakSites" <- peakSites
  processedFootprintData$"rawPeakFootprintMetrics" <- rawPeakFootprintMetrics
  processedFootprintData$"peakMotifSignal" <- peakMotifSignal
  processedFootprintData$"peakFlankSignal" <- peakFlankSignal
  processedFootprintData$"peakBackgroundSignal" <- peakBackgroundSignal
  processedFootprintData$"peak.log2Flank" <- peak.log2Flank
  processedFootprintData$"peak.log2Depth" <- peak.log2Depth
  ##
  processedFootprintData$"boundSites" <- boundSites
  processedFootprintData$"boundSitesMetrics" <- boundSitesMetrics
  processedFootprintData$"boundMotifSignal" <- boundMotifSignal
  processedFootprintData$"boundFlankSignal" <- boundFlankSignal
  processedFootprintData$"boundBackgroundSignal" <- boundBackgroundSignal
  processedFootprintData$"bound.log2Flank" <- bound.log2Flank
  processedFootprintData$"bound.log2Depth" <- bound.log2Depth
  ##
  processedFootprintData$"unboundSites" <- unboundSites
  processedFootprintData$"unboundSitesMetrics" <- unboundSitesMetrics
  processedFootprintData$"unboundMotifSignal" <- unboundMotifSignal
  processedFootprintData$"unboundFlankSignal" <- unboundFlankSignal
  processedFootprintData$"unboundBackgroundSignal" <- unboundBackgroundSignal
  processedFootprintData$"unbound.log2Flank" <- unbound.log2Flank
  processedFootprintData$"unbound.log2Depth" <- unbound.log2Depth
  
  ## Save the data
  save(processedFootprintData, file = dataOutPath)
  
  } # end if (numMotifs == 0)
} # end if (file.exists(dataOutPath) == TRUE)

}, # end try
error=function(cond){
  message(cond)
  return(NA)
},
finally={})
  
##
file.create(outPath)
cat("Finished!", "\n")