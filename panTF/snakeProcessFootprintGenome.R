
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
inputPath <- gsub("operations", "data", inputPath)
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
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
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
    cat("Setting genome sites", "\n")
    com <- paste0("genomeSites <- footprintData[['", motifNames[1], "']][['genomeSites']]")
    eval(parse(text = com))
    cat("Setting genome insertion matrix", "\n")
    com <- paste0("genomeInsertionMatrix <- footprintData[['", motifNames[1], "']][['insMatrix']]")
    eval(parse(text = com))
    cat("Setting footprint metrics", "\n")
    com <- paste0("rawGenomeFootprintMetrics <- footprintData[['", motifNames[1], "']][['rawFootprintMetrics']]")
    eval(parse(text = com))
    numGenomeSites <- length(genomeSites)
    cat("Number of genome sites", numGenomeSites, "\n")
    ## Split the sites into bound and unbound as determined by null model with bonferroni correction
    cat("Finding bound site overlaps", "\n")
    com <- paste0("boundSiteOverlaps <- findOverlaps(footprintData[['", motifNames[1], "']][['parseData']][['bfSites']], footprintData[['", motifNames[1], "']][['genomeSites']])")
    eval(parse(text = com))
    cat("Setting index for bound sites", "\n")
    boundSiteIndex <- boundSiteOverlaps@to
    cat("Subsetting bound sites", "\n")
    boundSites <- genomeSites[boundSiteIndex]
    boundSitesInsertionMatrix <- genomeInsertionMatrix[boundSiteIndex,]
    boundSitesMetrics <- rawGenomeFootprintMetrics[boundSiteIndex,]
    numBoundSites <- length(boundSites)
    cat("Subsetting unbound sites", "\n")
    unboundSites <- genomeSites[-boundSiteIndex]
    unboundSitesInsertionMatrix <- genomeInsertionMatrix[-boundSiteIndex,]
    unboundSitesMetrics <- rawGenomeFootprintMetrics[-boundSiteIndex,]
    numUnboundSites <- length(unboundSites)
    
    #### Calculate footprint characteristics on merged data ####
    ## Calculate the 10% trimmed mean of all insertions in the motif sites
    
    ## I am going to add in code here to convert these values to a per bp average
    ## Should move the code to parse script at some point
    com <- paste0("motifWidth <- footprintData[['", motifNames[1], "']][['motifWidth']]")
    eval(parse(text = com))
    ##
    genomeMotifSignal <- (mean(rawGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
    boundMotifSignal <- (mean(boundSitesMetrics[,3], trim = 0.10) / motifWidth)
    unboundMotifSignal <- (mean(unboundSitesMetrics[,3], trim = 0.10) / motifWidth)
    
    ## Calculate the mean of all insertions in the flank region
    genomeFlankSignal <- (mean(rawGenomeFootprintMetrics[,2], trim = 0.10) / 100)
    boundFlankSignal <- (mean(boundSitesMetrics[,2], trim = 0.10) / 100)
    unboundFlankSignal <- (mean(unboundSitesMetrics[,2], trim = 0.10) / 100)
    
    ## Calculate the mean of background insertions
    genomeBackgroundSignal <- (mean(rawGenomeFootprintMetrics[,1], trim = 0.10) / 100)
    boundBackgroundSignal <- (mean(boundSitesMetrics[,1], trim = 0.10) / 100)
    unboundBackgroundSignal <- (mean(unboundSitesMetrics[,1], trim = 0.10) / 100)
    
    ## Calculate flanking accessibility (log2 fold change between flank and background)
    genome.log2Flank <- log2(genomeFlankSignal / genomeBackgroundSignal)
    bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
    unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
    
    ## Calculate footprint depth (log2 fold change between flank and background)
    genome.log2Depth <- log2(genomeMotifSignal / genomeFlankSignal)
    bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
    unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
    
    #### TRANSFER DATA TO STORAGE OBJECT ####
    ## Initialize a new list object to store the processed data
    processedFootprintData <- list()
    ##
    processedFootprintData$"geneName" <- geneName
    processedFootprintData$"numMotifs" <- numMotifs
    processedFootprintData$"numGenomeSites" <- numGenomeSites
    processedFootprintData$"numBoundSites" <- numBoundSites
    processedFootprintData$"numUnboundSites" <- numUnboundSites
    ##
    processedFootprintData$"genomeSites" <- genomeSites
    processedFootprintData$"rawGenomeFootprintMetrics" <- rawGenomeFootprintMetrics
    processedFootprintData$"genomeMotifSignal" <- genomeMotifSignal
    processedFootprintData$"genomeFlankSignal" <- genomeFlankSignal
    processedFootprintData$"genomeBackgroundSignal" <- genomeBackgroundSignal
    processedFootprintData$"genome.log2Flank" <- genome.log2Flank
    processedFootprintData$"genome.log2Depth" <- genome.log2Depth
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
    
    
    #### CODE TESTING - PROMOTER/DISTAL groups ####
    ## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
    promoters <- promoters(txdb, upstream = 1000, downstream = 100)
    ## Trim the GRanges object to keep standard entries only
    scope <- paste0("chr", c(1:22, "X", "Y"))
    promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
    promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
    
    ## Subset based on the overlaps
    promoterOverlaps <- findOverlaps(promoters, genomeSites, ignore.strand = TRUE)
    promoterIdx <- unique(promoterOverlaps@to)
    ##
    promoterGenomeSites <- genomeSites[promoterIdx]
    distalGenomeSites <- genomeSites[-promoterIdx]
    ##
    promoterBoundOverlaps <- findOverlaps(promoterGenomeSites, boundSites)
    promoterUnboundOverlaps <- findOverlaps(promoterGenomeSites, unboundSites)
    promoterBoundIdx <- unique(promoterBoundOverlaps@from)
    promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
    ##
    promoterBoundSites <- genomeSites[promoterBoundIdx]
    promoterUnboundSites <- genomeSites[promoterUnboundIdx]
    ##
    distalBoundOverlaps <- findOverlaps(distalGenomeSites, boundSites)
    distalUnboundOverlaps <- findOverlaps(distalGenomeSites, unboundSites)
    distalBoundIdx <- unique(distalBoundOverlaps@from)
    distalUnboundIdx <- unique(distalUnboundOverlaps@from)
    ##
    distalBoundSites <- genomeSites[distalBoundIdx]
    distalUnboundSites <- genomeSites[distalUnboundIdx]
    
    ##### promoter genome sites
    promoterGenomeFootprintMetrics <- rawGenomeFootprintMetrics[promoterIdx,]
    promoterGenomeMotifSignal <- (mean(promoterGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
    promoterGenomeFlankSignal <- (mean(promoterGenomeFootprintMetrics[,2], trim = 0.10) / 100)
    promoterGenomeBackgroundSignal <- (mean(promoterGenomeFootprintMetrics[,1], trim = 0.10) / 100)
    promoterGenome.log2Flank <- log2(promoterGenomeFlankSignal / promoterGenomeBackgroundSignal)
    promoterGenome.log2Depth <- log2(promoterGenomeMotifSignal / promoterGenomeFlankSignal)
    
    ## distal genome Sites
    distalGenomeFootprintMetrics <- rawGenomeFootprintMetrics[-promoterIdx,]
    distalGenomeMotifSignal <- (mean(distalGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
    distalGenomeFlankSignal <- (mean(distalGenomeFootprintMetrics[,2], trim = 0.10) / 100)
    distalGenomeBackgroundSignal <- (mean(distalGenomeFootprintMetrics[,1], trim = 0.10) / 100)
    distalGenome.log2Flank <- log2(distalGenomeFlankSignal / distalGenomeBackgroundSignal)
    distalGenome.log2Depth <- log2(distalGenomeMotifSignal / distalGenomeFlankSignal)
    
    ## promoter Bound
    promoterBoundFootprintMetrics <- rawGenomeFootprintMetrics[promoterBoundIdx,]
    promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
    promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
    promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
    promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
    promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)
    
    ## promoter unbound
    promoterUnboundFootprintMetrics <- rawGenomeFootprintMetrics[promoterUnboundIdx,]
    promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
    promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
    promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
    promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
    promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)
    
    ## distal Bound
    distalBoundFootprintMetrics <- rawGenomeFootprintMetrics[distalBoundIdx,]
    distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
    distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
    distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
    distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
    distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)
    
    ## distal unbound
    distalUnboundFootprintMetrics <- rawGenomeFootprintMetrics[distalUnboundIdx,]
    distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
    distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
    distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
    distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
    distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)
    
    ## STORE THE DATA
    processedFootprintData$"promoterGenomeSites" <- promoterGenomeSites
    processedFootprintData$"promoterGenomeFootprintMetrics" <- promoterGenomeFootprintMetrics
    processedFootprintData$"promoterGenomeMotifSignal" <- promoterGenomeMotifSignal
    processedFootprintData$"promoterGenomeFlankSignal" <- promoterGenomeFlankSignal
    processedFootprintData$"promoterGenomeBackgroundSignal" <- promoterGenomeBackgroundSignal
    processedFootprintData$"promoterGenome.log2Flank" <- promoterGenome.log2Flank
    processedFootprintData$"promoterGenome.log2Depth" <- promoterGenome.log2Depth
    ##
    processedFootprintData$"distalGenomeSites" <- distalGenomeSites
    processedFootprintData$"distalGenomeFootprintMetrics" <- distalGenomeFootprintMetrics
    processedFootprintData$"distalGenomeMotifSignal" <- distalGenomeMotifSignal
    processedFootprintData$"distalGenomeFlankSignal" <- distalGenomeFlankSignal
    processedFootprintData$"distalGenomeBackgroundSignal" <- distalGenomeBackgroundSignal
    processedFootprintData$"distalGenome.log2Flank" <- distalGenome.log2Flank
    processedFootprintData$"distalGenome.log2Depth" <- distalGenome.log2Depth
    ##
    processedFootprintData$"promoterBoundSites" <- promoterBoundSites
    processedFootprintData$"promoterBoundFootprintMetrics" <- promoterBoundFootprintMetrics
    processedFootprintData$"promoterBoundMotifSignal" <- promoterBoundMotifSignal
    processedFootprintData$"promoterBoundFlankSignal" <- promoterBoundFlankSignal
    processedFootprintData$"promoterBoundBackgroundSignal" <- promoterBoundBackgroundSignal
    processedFootprintData$"promoterBound.log2Flank" <- promoterBound.log2Flank
    processedFootprintData$"promoterBound.log2Depth" <- promoterBound.log2Depth
    ##
    processedFootprintData$"promoterUnboundSites" <- promoterUnboundSites
    processedFootprintData$"promoterUnboundFootprintMetrics" <- promoterUnboundFootprintMetrics
    processedFootprintData$"promoterUnboundMotifSignal" <- promoterUnboundMotifSignal
    processedFootprintData$"promoterUnboundFlankSignal" <- promoterUnboundFlankSignal
    processedFootprintData$"promoterUnboundBackgroundSignal" <- promoterUnboundBackgroundSignal
    processedFootprintData$"promoterUnbound.log2Flank" <- promoterUnbound.log2Flank
    processedFootprintData$"promoterUnbound.log2Depth" <- promoterUnbound.log2Depth
    ##
    processedFootprintData$"distalBoundSites" <- distalBoundSites
    processedFootprintData$"distalBoundFootprintMetrics" <- distalBoundFootprintMetrics
    processedFootprintData$"distalBoundMotifSignal" <- distalBoundMotifSignal
    processedFootprintData$"distalBoundFlankSignal" <- distalBoundFlankSignal
    processedFootprintData$"distalBoundBackgroundSignal" <- distalBoundBackgroundSignal
    processedFootprintData$"distalBound.log2Flank" <- distalBound.log2Flank
    processedFootprintData$"distalBound.log2Depth" <- distalBound.log2Depth
    ##
    processedFootprintData$"distalUnboundSites" <- distalUnboundSites
    processedFootprintData$"distalUnboundFootprintMetrics" <- distalUnboundFootprintMetrics
    processedFootprintData$"distalUnboundMotifSignal" <- distalUnboundMotifSignal
    processedFootprintData$"distalUnboundFlankSignal" <- distalUnboundFlankSignal
    processedFootprintData$"distalUnboundBackgroundSignal" <- distalUnboundBackgroundSignal
    processedFootprintData$"distalUnbound.log2Flank" <- distalUnbound.log2Flank
    processedFootprintData$"distalUnbound.log2Depth" <- distalUnbound.log2Depth
 
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
      cat("Pulling genome sites", "\n")
      com <- paste0("genomeSites", z, " <- footprintData[['", motifNames[z], "']][['genomeSites']]")
      eval(parse(text = com))
      cat("Pulling insertion matrix", "\n")
      com <- paste0("genomeInsertionMatrix", z, " <- footprintData[['", motifNames[z], "']][['insMatrix']]")
      eval(parse(text = com))
      cat("Pulling footprint metrics", "\n")
      com <- paste0("rawGenomeFootprintMetrics", z, " <- footprintData[['", motifNames[z], "']][['rawFootprintMetrics']]")
      eval(parse(text = com))
      ## Perform the bound/unbound split
      cat("Finding bound overlaps", "\n")
      com <- paste0("boundSiteOverlaps", z, " <- findOverlaps(footprintData[['", motifNames[z], "']][['parseData']][['bfSites']], footprintData[['", motifNames[z], "']][['genomeSites']])")
      eval(parse(text = com))
      cat("Setting bound index", "\n")
      com <- paste0("boundSiteIndex", z, " <- boundSiteOverlaps", z, "@to")
      eval(parse(text = com))
      ## Bound sites
      cat("Setting bound site data", "\n")
      com <- paste0("boundSites", z, " <- genomeSites", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("boundSitesInsertionMatrix", z, " <- genomeInsertionMatrix", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("boundSitesMetrics", z, " <- rawGenomeFootprintMetrics", z, "[boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("numBoundSites", z, " <- length(boundSites", z, ")")
      eval(parse(text = com))
      ## Unbound sites
      cat("Setting unbound site data", "\n")
      com <- paste0("unboundSites", z, " <- genomeSites", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("unboundSitesInsertionMatrix", z, " <- genomeInsertionMatrix", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("unboundSitesMetrics", z, " <- rawGenomeFootprintMetrics", z, "[-boundSiteIndex", z, ",]")
      eval(parse(text = com))
      com <- paste0("numUnboundSites", z, " <- length(unboundSites", z, ")")
      eval(parse(text = com))
    } # end for (z in 1:numMotif)
    
    #### FIX THIS CODE LATER ####
    #### ADJUSTING VALUES TO PER BP ####
    for (m in 1:numMotifs){
      
      com <- paste0("motifWidth <- footprintData[['", motifNames[m], "']][['motifWidth']]")
      eval(parse(text = com))
      
      com <- paste0("boundSitesMetrics", m, "[,1] <- boundSitesMetrics", m, "[,1] / 100")
      eval(parse(text = com))
      com <- paste0("boundSitesMetrics", m, "[,2] <- boundSitesMetrics", m, "[,2] / 100")
      eval(parse(text = com))
      com <- paste0("boundSitesMetrics", m, "[,3] <- boundSitesMetrics", m, "[,3] / motifWidth")
      eval(parse(text = com))
      
      com <- paste0("rawGenomeFootprintMetrics", m, "[,1] <- rawGenomeFootprintMetrics", m, "[,1] / 100")
      eval(parse(text = com))
      com <- paste0("rawGenomeFootprintMetrics", m, "[,2] <- rawGenomeFootprintMetrics", m, "[,2] / 100")
      eval(parse(text = com))
      com <- paste0("rawGenomeFootprintMetrics", m, "[,3] <- rawGenomeFootprintMetrics", m, "[,3] / motifWidth")
      eval(parse(text = com))
      
      com <- paste0("unboundSitesMetrics", m, "[,1] <- unboundSitesMetrics", m, "[,1] / 100")
      eval(parse(text = com))
      com <- paste0("unboundSitesMetrics", m, "[,2] <- unboundSitesMetrics", m, "[,2] / 100")
      eval(parse(text = com))
      com <- paste0("unboundSitesMetrics", m, "[,3] <- unboundSitesMetrics", m, "[,3] / motifWidth")
      eval(parse(text = com))
    }
    
    #### Perform the merging and deduplication ####
    for (b in 2:numMotifs){
      ## Find overlaps and generate selection indices
      com <- paste0("tempOverlapsGenome <- findOverlaps(genomeSites1, genomeSites", b,")")
      eval(parse(text = com))
      com <- paste0("tempOverlapsBound <- findOverlaps(boundSites1, boundSites", b,")")
      eval(parse(text = com))
      com <- paste0("tempOverlapsUnbound <- findOverlaps(unboundSites1, unboundSites", b,")")
      eval(parse(text = com))
      
      #### MERGE THE GENOME SITES (all sites) ####
      ## If no overlaps are present, can just directly merge the two Granges
      ## Otherwise, if some overlaps are present, merge the Granges, but omit overlapping sites from second group
      if (length(tempOverlapsGenome@from) == 0){
        com <- paste0("genomeSites1 <- c(genomeSites1, genomeSites", b, ")")
        eval(parse(text = com))
        com <- paste0("rawGenomeFootprintMetrics1 <- rbind(rawGenomeFootprintMetrics1, rawGenomeFootprintMetrics", b, ")")
        eval(parse(text = com))
      } else {
        mergeIdx <- tempOverlapsGenome@to
        com <- paste0("genomeSites1 <- c(genomeSites1, genomeSites", b, "[-mergeIdx])")
        eval(parse(text = com))
        com <- paste0("rawGenomeFootprintMetrics1 <- rbind(rawGenomeFootprintMetrics1, rawGenomeFootprintMetrics", b, "[-mergeIdx,])")
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
    genomeSites <- genomeSites1
    rawGenomeFootprintMetrics <- rawGenomeFootprintMetrics1
    numGenomeSites <- length(genomeSites)
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
  genomeMotifSignal <- mean(rawGenomeFootprintMetrics[,3], trim = 0.10)
  boundMotifSignal <- mean(boundSitesMetrics[,3], trim = 0.10)
  unboundMotifSignal <- mean(unboundSitesMetrics[,3], trim = 0.10)
  
  ## Calculate the mean of all insertions in the flank region
  genomeFlankSignal <- mean(rawGenomeFootprintMetrics[,2], trim = 0.10)
  boundFlankSignal <- mean(boundSitesMetrics[,2], trim = 0.10)
  unboundFlankSignal <- mean(unboundSitesMetrics[,2], trim = 0.10)
  
  ## Calculate the mean of background insertions
  genomeBackgroundSignal <- mean(rawGenomeFootprintMetrics[,1], trim = 0.10)
  boundBackgroundSignal <- mean(boundSitesMetrics[,1], trim = 0.10)
  unboundBackgroundSignal <- mean(unboundSitesMetrics[,1], trim = 0.10)
  
  ## Calculate flanking accessibility (log2 fold change between flank and background)
  genome.log2Flank <- log2(genomeFlankSignal / genomeBackgroundSignal)
  bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
  unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
  
  ## Calculate footprint depth (log2 fold change between flank and background)
  genome.log2Depth <- log2(genomeMotifSignal / genomeFlankSignal)
  bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
  unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
  
  #### TRANSFER DATA TO STORAGE OBJECT ####
  ## Initialize a new list object to store the processed data
  processedFootprintData <- list()
  ##
  processedFootprintData$"geneName" <- geneName
  processedFootprintData$"numMotifs" <- numMotifs
  processedFootprintData$"numGenomeSites" <- numGenomeSites
  processedFootprintData$"numBoundSites" <- numBoundSites
  processedFootprintData$"numUnboundSites" <- numUnboundSites
  ##
  processedFootprintData$"genomeSites" <- genomeSites
  processedFootprintData$"rawGenomeFootprintMetrics" <- rawGenomeFootprintMetrics
  processedFootprintData$"genomeMotifSignal" <- genomeMotifSignal
  processedFootprintData$"genomeFlankSignal" <- genomeFlankSignal
  processedFootprintData$"genomeBackgroundSignal" <- genomeBackgroundSignal
  processedFootprintData$"genome.log2Flank" <- genome.log2Flank
  processedFootprintData$"genome.log2Depth" <- genome.log2Depth
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
  
  #### CODE TESTING - PROMOTER/DISTAL groups ####
  ## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
  promoters <- promoters(txdb, upstream = 1000, downstream = 100)
  ## Trim the GRanges object to keep standard entries only
  scope <- paste0("chr", c(1:22, "X", "Y"))
  promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
  promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
  
  ## Subset based on the overlaps
  promoterOverlaps <- findOverlaps(promoters, genomeSites, ignore.strand = TRUE)
  promoterIdx <- unique(promoterOverlaps@to)
  ##
  promoterGenomeSites <- genomeSites[promoterIdx]
  distalGenomeSites <- genomeSites[-promoterIdx]
  ##
  promoterBoundOverlaps <- findOverlaps(promoterGenomeSites, boundSites)
  promoterUnboundOverlaps <- findOverlaps(promoterGenomeSites, unboundSites)
  promoterBoundIdx <- unique(promoterBoundOverlaps@from)
  promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
  ##
  promoterBoundSites <- genomeSites[promoterBoundIdx]
  promoterUnboundSites <- genomeSites[promoterUnboundIdx]
  ##
  distalBoundOverlaps <- findOverlaps(distalGenomeSites, boundSites)
  distalUnboundOverlaps <- findOverlaps(distalGenomeSites, unboundSites)
  distalBoundIdx <- unique(distalBoundOverlaps@from)
  distalUnboundIdx <- unique(distalUnboundOverlaps@from)
  ##
  distalBoundSites <- genomeSites[distalBoundIdx]
  distalUnboundSites <- genomeSites[distalUnboundIdx]
  
  ##### promoter genome sites
  promoterGenomeFootprintMetrics <- rawGenomeFootprintMetrics[promoterIdx,]
  promoterGenomeMotifSignal <- (mean(promoterGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterGenomeFlankSignal <- (mean(promoterGenomeFootprintMetrics[,2], trim = 0.10) / 100)
  promoterGenomeBackgroundSignal <- (mean(promoterGenomeFootprintMetrics[,1], trim = 0.10) / 100)
  promoterGenome.log2Flank <- log2(promoterGenomeFlankSignal / promoterGenomeBackgroundSignal)
  promoterGenome.log2Depth <- log2(promoterGenomeMotifSignal / promoterGenomeFlankSignal)
  
  ## distal genome Sites
  distalGenomeFootprintMetrics <- rawGenomeFootprintMetrics[-promoterIdx,]
  distalGenomeMotifSignal <- (mean(distalGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalGenomeFlankSignal <- (mean(distalGenomeFootprintMetrics[,2], trim = 0.10) / 100)
  distalGenomeBackgroundSignal <- (mean(distalGenomeFootprintMetrics[,1], trim = 0.10) / 100)
  distalGenome.log2Flank <- log2(distalGenomeFlankSignal / distalGenomeBackgroundSignal)
  distalGenome.log2Depth <- log2(distalGenomeMotifSignal / distalGenomeFlankSignal)
  
  ## promoter Bound
  promoterBoundFootprintMetrics <- rawGenomeFootprintMetrics[promoterBoundIdx,]
  promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
  promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
  promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
  promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)
  
  ## promoter unbound
  promoterUnboundFootprintMetrics <- rawGenomeFootprintMetrics[promoterUnboundIdx,]
  promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
  promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
  promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
  promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)
  
  ## distal Bound
  distalBoundFootprintMetrics <- rawGenomeFootprintMetrics[distalBoundIdx,]
  distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
  distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
  distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
  distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)
  
  ## distal unbound
  distalUnboundFootprintMetrics <- rawGenomeFootprintMetrics[distalUnboundIdx,]
  distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
  distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
  distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
  distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)
  
  ## STORE THE DATA
  processedFootprintData$"promoterGenomeSites" <- promoterGenomeSites
  processedFootprintData$"promoterGenomeFootprintMetrics" <- promoterGenomeFootprintMetrics
  processedFootprintData$"promoterGenomeMotifSignal" <- promoterGenomeMotifSignal
  processedFootprintData$"promoterGenomeFlankSignal" <- promoterGenomeFlankSignal
  processedFootprintData$"promoterGenomeBackgroundSignal" <- promoterGenomeBackgroundSignal
  processedFootprintData$"promoterGenome.log2Flank" <- promoterGenome.log2Flank
  processedFootprintData$"promoterGenome.log2Depth" <- promoterGenome.log2Depth
  ##
  processedFootprintData$"distalGenomeSites" <- distalGenomeSites
  processedFootprintData$"distalGenomeFootprintMetrics" <- distalGenomeFootprintMetrics
  processedFootprintData$"distalGenomeMotifSignal" <- distalGenomeMotifSignal
  processedFootprintData$"distalGenomeFlankSignal" <- distalGenomeFlankSignal
  processedFootprintData$"distalGenomeBackgroundSignal" <- distalGenomeBackgroundSignal
  processedFootprintData$"distalGenome.log2Flank" <- distalGenome.log2Flank
  processedFootprintData$"distalGenome.log2Depth" <- distalGenome.log2Depth
  ##
  processedFootprintData$"promoterBoundSites" <- promoterBoundSites
  processedFootprintData$"promoterBoundFootprintMetrics" <- promoterBoundFootprintMetrics
  processedFootprintData$"promoterBoundMotifSignal" <- promoterBoundMotifSignal
  processedFootprintData$"promoterBoundFlankSignal" <- promoterBoundFlankSignal
  processedFootprintData$"promoterBoundBackgroundSignal" <- promoterBoundBackgroundSignal
  processedFootprintData$"promoterBound.log2Flank" <- promoterBound.log2Flank
  processedFootprintData$"promoterBound.log2Depth" <- promoterBound.log2Depth
  ##
  processedFootprintData$"promoterUnboundSites" <- promoterUnboundSites
  processedFootprintData$"promoterUnboundFootprintMetrics" <- promoterUnboundFootprintMetrics
  processedFootprintData$"promoterUnboundMotifSignal" <- promoterUnboundMotifSignal
  processedFootprintData$"promoterUnboundFlankSignal" <- promoterUnboundFlankSignal
  processedFootprintData$"promoterUnboundBackgroundSignal" <- promoterUnboundBackgroundSignal
  processedFootprintData$"promoterUnbound.log2Flank" <- promoterUnbound.log2Flank
  processedFootprintData$"promoterUnbound.log2Depth" <- promoterUnbound.log2Depth
  ##
  processedFootprintData$"distalBoundSites" <- distalBoundSites
  processedFootprintData$"distalBoundFootprintMetrics" <- distalBoundFootprintMetrics
  processedFootprintData$"distalBoundMotifSignal" <- distalBoundMotifSignal
  processedFootprintData$"distalBoundFlankSignal" <- distalBoundFlankSignal
  processedFootprintData$"distalBoundBackgroundSignal" <- distalBoundBackgroundSignal
  processedFootprintData$"distalBound.log2Flank" <- distalBound.log2Flank
  processedFootprintData$"distalBound.log2Depth" <- distalBound.log2Depth
  ##
  processedFootprintData$"distalUnboundSites" <- distalUnboundSites
  processedFootprintData$"distalUnboundFootprintMetrics" <- distalUnboundFootprintMetrics
  processedFootprintData$"distalUnboundMotifSignal" <- distalUnboundMotifSignal
  processedFootprintData$"distalUnboundFlankSignal" <- distalUnboundFlankSignal
  processedFootprintData$"distalUnboundBackgroundSignal" <- distalUnboundBackgroundSignal
  processedFootprintData$"distalUnbound.log2Flank" <- distalUnbound.log2Flank
  processedFootprintData$"distalUnbound.log2Depth" <- distalUnbound.log2Depth
  
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