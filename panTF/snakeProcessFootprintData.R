

##########
rawInsProb
rawTotalSignal





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

## I am going to add in code here to convert these values to a per bp average
## Should move the code to parse script at some point
com <- paste0("motifWidth <- footprintData[['", motifNames[1], "']][['motifWidth']]")
eval(parse(text = com))
##
peakMotifSignal <- (mean(rawPeakFootprintMetrics[,3], trim = 0.10) / motifWidth)
boundMotifSignal <- (mean(boundSitesMetrics[,3], trim = 0.10) / motifWidth)
unboundMotifSignal <- (mean(unboundSitesMetrics[,3], trim = 0.10) / motifWidth)

## Calculate the mean of all insertions in the flank region
peakFlankSignal <- (mean(rawPeakFootprintMetrics[,2], trim = 0.10) / 100)
boundFlankSignal <- (mean(boundSitesMetrics[,2], trim = 0.10) / 100)
unboundFlankSignal <- (mean(unboundSitesMetrics[,2], trim = 0.10) / 100)

## Calculate the mean of background insertions
peakBackgroundSignal <- (mean(rawPeakFootprintMetrics[,1], trim = 0.10) / 100)
boundBackgroundSignal <- (mean(boundSitesMetrics[,1], trim = 0.10) / 100)
unboundBackgroundSignal <- (mean(unboundSitesMetrics[,1], trim = 0.10) / 100)

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
processedFootprintData$"geneName" <- geneName
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


#### CODE TESTING - PROMOTER/DISTAL groups ####
## Pull promoters from txdb, define promoter region as -1000/+500 in accordance with TCGA paper
promoters <- promoters(txdb, upstream = 1000, downstream = 100)
## Trim the GRanges object to keep standard entries only
scope <- paste0("chr", c(1:22, "X", "Y"))
promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")

## Subset based on the overlaps
promoterOverlaps <- findOverlaps(promoters, peakSites, ignore.strand = TRUE)
promoterIdx <- unique(promoterOverlaps@to)
##
promoterPeakSites <- peakSites[promoterIdx]
distalPeakSites <- peakSites[-promoterIdx]
##
promoterBoundOverlaps <- findOverlaps(promoterPeakSites, boundSites)
promoterUnboundOverlaps <- findOverlaps(promoterPeakSites, unboundSites)
promoterBoundIdx <- unique(promoterBoundOverlaps@from)
promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
##
promoterBoundSites <- peakSites[promoterBoundIdx]
promoterUnboundSites <- peakSites[promoterUnboundIdx]
##
distalBoundOverlaps <- findOverlaps(distalPeakSites, boundSites)
distalUnboundOverlaps <- findOverlaps(distalPeakSites, unboundSites)
distalBoundIdx <- unique(distalBoundOverlaps@from)
distalUnboundIdx <- unique(distalUnboundOverlaps@from)
##
distalBoundSites <- peakSites[distalBoundIdx]
distalUnboundSites <- peakSites[distalUnboundIdx]

##### promoter peak sites
promoterPeakFootprintMetrics <- rawPeakFootprintMetrics[promoterIdx,]
promoterPeakMotifSignal <- (mean(promoterPeakFootprintMetrics[,3], trim = 0.10) / motifWidth)
promoterPeakFlankSignal <- (mean(promoterPeakFootprintMetrics[,2], trim = 0.10) / 100)
promoterPeakBackgroundSignal <- (mean(promoterPeakFootprintMetrics[,1], trim = 0.10) / 100)
promoterPeak.log2Flank <- log2(promoterPeakFlankSignal / promoterPeakBackgroundSignal)
promoterPeak.log2Depth <- log2(promoterPeakMotifSignal / promoterPeakFlankSignal)

## distalPeakSites
distalPeakFootprintMetrics <- rawPeakFootprintMetrics[-promoterIdx,]
distalPeakMotifSignal <- (mean(distalPeakFootprintMetrics[,3], trim = 0.10) / motifWidth)
distalPeakFlankSignal <- (mean(distalPeakFootprintMetrics[,2], trim = 0.10) / 100)
distalPeakBackgroundSignal <- (mean(distalPeakFootprintMetrics[,1], trim = 0.10) / 100)
distalPeak.log2Flank <- log2(distalPeakFlankSignal / distalPeakBackgroundSignal)
distalPeak.log2Depth <- log2(distalPeakMotifSignal / distalPeakFlankSignal)

## promoter Bound
promoterBoundFootprintMetrics <- rawPeakFootprintMetrics[promoterBoundIdx,]
promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)

## promoter unbound
promoterUnboundFootprintMetrics <- rawPeakFootprintMetrics[promoterUnboundIdx,]
promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)

## distal Bound
distalBoundFootprintMetrics <- rawPeakFootprintMetrics[distalBoundIdx,]
distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)

## distal unbound
distalUnboundFootprintMetrics <- rawPeakFootprintMetrics[distalUnboundIdx,]
distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)

## STORE THE DATA
processedFootprintData$"promoterPeakSites" <- promoterPeakSites
processedFootprintData$"promoterPeakFootprintMetrics" <- promoterPeakFootprintMetrics
processedFootprintData$"promoterPeakMotifSignal" <- promoterPeakMotifSignal
processedFootprintData$"promoterPeakFlankSignal" <- promoterPeakFlankSignal
processedFootprintData$"promoterPeakBackgroundSignal" <- promoterPeakBackgroundSignal
processedFootprintData$"promoterPeak.log2Flank" <- promoterPeak.log2Flank
processedFootprintData$"promoterPeak.log2Depth" <- promoterPeak.log2Depth
##
processedFootprintData$"distalpeakSites" <- distalPeakSites
processedFootprintData$"distalPeakFootprintMetrics" <- distalPeakFootprintMetrics
processedFootprintData$"distalPeakMotifSignal" <- distalPeakMotifSignal
processedFootprintData$"distalPeakFlankSignal" <- distalPeakFlankSignal
processedFootprintData$"distalPeakBackgroundSignal" <- distalPeakBackgroundSignal
processedFootprintData$"distalPeak.log2Flank" <- distalPeak.log2Flank
processedFootprintData$"distalPeak.log2Depth" <- distalPeak.log2Depth
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
    
    com <- paste0("rawPeakFootprintMetrics", m, "[,1] <- rawPeakFootprintMetrics", m, "[,1] / 100")
    eval(parse(text = com))
    com <- paste0("rawPeakFootprintMetrics", m, "[,2] <- rawPeakFootprintMetrics", m, "[,2] / 100")
    eval(parse(text = com))
    com <- paste0("rawPeakFootprintMetrics", m, "[,3] <- rawPeakFootprintMetrics", m, "[,3] / motifWidth")
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
  processedFootprintData$"geneName" <- geneName
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
  
  #### CODE TESTING - PROMOTER/DISTAL groups ####
  ## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
  promoters <- promoters(txdb, upstream = 1000, downstream = 100)
  ## Trim the GRanges object to keep standard entries only
  scope <- paste0("chr", c(1:22, "X", "Y"))
  promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
  promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
  
  ## Subset based on the overlaps
  promoterOverlaps <- findOverlaps(promoters, peakSites, ignore.strand = TRUE)
  promoterIdx <- unique(promoterOverlaps@to)
  ##
  promoterPeakSites <- peakSites[promoterIdx]
  distalPeakSites <- peakSites[-promoterIdx]
  ##
  promoterBoundOverlaps <- findOverlaps(promoterPeakSites, boundSites)
  promoterUnboundOverlaps <- findOverlaps(promoterPeakSites, unboundSites)
  promoterBoundIdx <- unique(promoterBoundOverlaps@from)
  promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
  ##
  promoterBoundSites <- peakSites[promoterBoundIdx]
  promoterUnboundSites <- peakSites[promoterUnboundIdx]
  ##
  distalBoundOverlaps <- findOverlaps(distalPeakSites, boundSites)
  distalUnboundOverlaps <- findOverlaps(distalPeakSites, unboundSites)
  distalBoundIdx <- unique(distalBoundOverlaps@from)
  distalUnboundIdx <- unique(distalUnboundOverlaps@from)
  ##
  distalBoundSites <- peakSites[distalBoundIdx]
  distalUnboundSites <- peakSites[distalUnboundIdx]
  
  ##### promoter peak sites
  promoterPeakFootprintMetrics <- rawPeakFootprintMetrics[promoterIdx,]
  promoterPeakMotifSignal <- (mean(promoterPeakFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterPeakFlankSignal <- (mean(promoterPeakFootprintMetrics[,2], trim = 0.10) / 100)
  promoterPeakBackgroundSignal <- (mean(promoterPeakFootprintMetrics[,1], trim = 0.10) / 100)
  promoterPeak.log2Flank <- log2(promoterPeakFlankSignal / promoterPeakBackgroundSignal)
  promoterPeak.log2Depth <- log2(promoterPeakMotifSignal / promoterPeakFlankSignal)
  
  ## distalPeakSites
  distalPeakFootprintMetrics <- rawPeakFootprintMetrics[-promoterIdx,]
  distalPeakMotifSignal <- (mean(distalPeakFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalPeakFlankSignal <- (mean(distalPeakFootprintMetrics[,2], trim = 0.10) / 100)
  distalPeakBackgroundSignal <- (mean(distalPeakFootprintMetrics[,1], trim = 0.10) / 100)
  distalPeak.log2Flank <- log2(distalPeakFlankSignal / distalPeakBackgroundSignal)
  distalPeak.log2Depth <- log2(distalPeakMotifSignal / distalPeakFlankSignal)
  
  ## promoter Bound
  promoterBoundFootprintMetrics <- rawPeakFootprintMetrics[promoterBoundIdx,]
  promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
  promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
  promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
  promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)
  
  ## promoter unbound
  promoterUnboundFootprintMetrics <- rawPeakFootprintMetrics[promoterUnboundIdx,]
  promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
  promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
  promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
  promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)
  
  ## distal Bound
  distalBoundFootprintMetrics <- rawPeakFootprintMetrics[distalBoundIdx,]
  distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
  distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
  distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
  distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)
  
  ## distal unbound
  distalUnboundFootprintMetrics <- rawPeakFootprintMetrics[distalUnboundIdx,]
  distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
  distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
  distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
  distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
  distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)
  
  ## STORE THE DATA
  processedFootprintData$"promoterPeakSites" <- promoterPeakSites
  processedFootprintData$"promoterPeakFootprintMetrics" <- promoterPeakFootprintMetrics
  processedFootprintData$"promoterPeakMotifSignal" <- promoterPeakMotifSignal
  processedFootprintData$"promoterPeakFlankSignal" <- promoterPeakFlankSignal
  processedFootprintData$"promoterPeakBackgroundSignal" <- promoterPeakBackgroundSignal
  processedFootprintData$"promoterPeak.log2Flank" <- promoterPeak.log2Flank
  processedFootprintData$"promoterPeak.log2Depth" <- promoterPeak.log2Depth
  ##
  processedFootprintData$"distalpeakSites" <- distalPeakSites
  processedFootprintData$"distalPeakFootprintMetrics" <- distalPeakFootprintMetrics
  processedFootprintData$"distalPeakMotifSignal" <- distalPeakMotifSignal
  processedFootprintData$"distalPeakFlankSignal" <- distalPeakFlankSignal
  processedFootprintData$"distalPeakBackgroundSignal" <- distalPeakBackgroundSignal
  processedFootprintData$"distalPeak.log2Flank" <- distalPeak.log2Flank
  processedFootprintData$"distalPeak.log2Depth" <- distalPeak.log2Depth
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
################################################################################################################

## Disable scientific notation in variables
options(scipen = 999)
## suppress warnings globally here, as they will disrupt the tryCatch block
## will need to improve this code at some point
options(warn = -1)

#### Set snakemake variables
footprintDataPath <- snakemake@input[[1]]
sampleTotalReadsPath <- snakemake@input[[2]]
peaksPath <- snakemake@input[[3]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

#### Set hg38 number of bases (haploid)
hg38TotalBP <- 3272116950

#### Set the output filepath for the Rdata object and perform a filecheck
cat("Processing footprint data for", geneName, "from sample", sampleName, "\n")
dataOutPath <- gsub("operations/parse", "data/parsed", outPath)
dataOutPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", dataOutPath, perl = TRUE)
cat("Output path for parsed data:", dataOutPath, "\n")

#### Perform a filecheck 
if (file.exists(dataOutPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  #### Load libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(stats4))
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(parallel))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(genomation))
  suppressMessages(library(seqLogo))
  suppressMessages(library(ChIPpeakAnno))
  suppressMessages(library(rlist))
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  #### Load hg38 annotations
  cat("Loading hg38 annotations from TxDb", "\n")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
  
  #### Build functions
  cat("Building functions...", "\n")
  generateNullFP <- function(iterations, inputSignal, analysisWidth, motifWidth){
    # This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
    # Consider the total signal (number of insertions) at each specific ~200 bp locus
    # Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
    # Use published or experimentally derived models of Tn5 sequence specific insertion bias
    # For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
    # Generate the null model graph by weighted random residstribution of the total observed signal at that site
    # Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
    # These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
    # iterations = number of iterations
    # inputSignals = unique values for total signal
    # analysisWidth = total bp in region of interest (flank + background + motif)
    # motifWidth = motif width
    
    # declare vector of size n to store average motif signal values
    averages <- c()
    
    # generate the null models and calculate motif averages
    for (a in 1:iterations){
      
      # declare the null vector
      null <- c(1:(analysisWidth))
      
      # randomly distribute the total signal
      # size = the number of values to distribute
      # prob = probability of each site
      # length = length of the generated vector
      null <- c(as.vector(rmultinom(1, size=inputSignal, prob=rep(1, length(null)))))
      
      ## Calculate the mean signal in motif region
      motifStart <- ((analysisWidth - motifWidth)/2)
      motifEnd <- (motifStart + motifWidth)
      motifAvg <- (sum(null[motifStart:motifEnd])) / motifWidth
      
      ## Store the average values
      averages[a] <- motifAvg
      
    } # end for (a in 1:n)
    return(averages)
  } # end generateNullFP function
  
  #### Load the raw footprintData file
  footprintDataPath <- gsub("operations", "data", footprintDataPath)
  footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
  cat("Loading raw footprintData file from path:", footprintDataPath, "\n")
  load(footprintDataPath)
  
  #### To avoid errors, clear the list of any empty sub-lists first
  ## Remove motifs with 0 binding sites
  cat("Removing motifs that matched 0 genomic loci", "\n")
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  numMotif <- length(footprintData)
  cat("Found data for", numMotif, "motifs", "\n")
  
  #### If no raw data is found, skip this footprint
  if (numMotif == 0){
    cat(numMotif, "No data to analyze. Skipping", "\n")
  } else {
    
    #### Get the total number of reads in the sample
    cat("Loading total sample reads from:", sampleTotalReadsPath, "\n")
    load(sampleTotalReadsPath)
    cat("Found", sampleTotalReads, "total reads in current sample", "\n")
    
    #### Load the peaks data for current sample
    cat("Loading sample accesibility peak data from:", peaksPath, "\n")
    samplePeaks <- readBed(peaksPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
    samplePeaks <- keepStandardChromosomes(samplePeaks, pruning.mode="coarse")
    samplePeaks <- trim(samplePeaks)
    
    #### Loop over all motifs, perform the processing operations
    cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
    for (a in 1:numMotif){
      
      #### Pull the data from the raw footprint analysis
      cat("Processing footprint data for motif", a, "\n")
      com <- paste0("tempData <- footprintData$motif", a)
      eval(parse(text = com))
      ##
      PWM <- tempData[["PWM"]]
      motifWidth <- tempData[["motifWidth"]]
      allSites <- tempData[["allSites"]]
      extSites <- tempData[["extSites"]]
      insMatrix <- tempData[["insMatrix"]]
      libSize <- tempData[["libsize"]]
      coverageSize <- tempData[["coverageSize"]]
      libFactor <- tempData[["libFactor"]]
      ## Note that because trimming to standard xsomes is done, will need to set
      ## the total number of sites to the row number of the insertion matrix here
      numSites <- length(insMatrix[,1])
      
      #### Calculate basic statistics for each site with raw data
      rawSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 10)
      colnames(rawSiteBasicStats) <- c("Site index", "Total signal", "Total signal per bp", "Motif signal per bp",
                                    "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                    "Flank / Background", "Motif / Flank", "Motif / Wide Flank") 
      
      #### Populate the basic stats matrix
      for (b in 1:numSites){
        rawSiteBasicStats[b,1] <- b # Site index
        rawSiteBasicStats[b,2] <- sum(insMatrix[b,]) # Total signal
        rawSiteBasicStats[b,3] <- rawSiteBasicStats[b,2] / (500 + motifWidth) # Total signal per bp
        rawSiteBasicStats[b,4] <- sum(insMatrix[b,(250:(250 + motifWidth))] / motifWidth) # Motif signal per bp
        rawSiteBasicStats[b,5] <- (sum(insMatrix[b,200:250]) + sum(insMatrix[b,(250 + motifWidth):(300 + motifWidth)])) / 100 # Flank signal per bp
        rawSiteBasicStats[b,6] <- (sum(insMatrix[b,1:50]) + sum(insMatrix[b,(500 + motifWidth-50):(500 + motifWidth)])) / 100 # Background signal per bp
        rawSiteBasicStats[b,7] <- (sum(insMatrix[b,1:250]) + sum(insMatrix[b,(250 + motifWidth):(500 + motifWidth)])) / 500 # Wide flank signal per bp
        rawSiteBasicStats[b,8] <- rawSiteBasicStats[b,5] / rawSiteBasicStats[b,6] # Flank / background
        rawSiteBasicStats[b,9] <- rawSiteBasicStats[b,4] / rawSiteBasicStats[b,5] # Motif / flank
        rawSiteBasicStats[b,10] <- rawSiteBasicStats[b,4] / rawSiteBasicStats[b,7] # Motif / wide flank
      } # end for (b in 1:numSites)
      
      #### Generate the insertion site probability vector for raw data
      rawInsProb <- c()
      for (c in 1:(500 + motifWidth)){
        rawInsProb[c] <- sum(insMatrix[,c])
      } # end for (c in 1:(500 + motifWidth))
      rawTotalSignal<- sum(rawInsProb)
      rawInsProb <- rawInsProb / rawTotalSignal
        
      
      
      
      ## NEED TO CHECK MATH ON THIS ##
      #### Normalize values in insertion matrix with z-scores
      insStandardDeviation <- sd(insMatrix)
      insMean <- mean(insMatrix)
      ##
      zscoreInsMatrix <- ((insMatrix - insMean) / insStandardDeviation)
      
      ## Calculate basic statistics for z-score normalized matrix
      zscoreBasicStats <- matrix(data = NA, nrow = numSites, ncol = 10)
      colnames(zscoreBasicStats) <- c("Site index", "Total signal", "Total signal per bp", "Motif signal per bp",
                                      "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                      "Flank / Background", "Motif / Flank", "Motif / Wide Flank") 
      
      ## Populate the zscore stats matrix
      for (b in 1:numSites){
        # Site index
        zscoreBasicStats[b,1] <- b
        # Total signal
        zscoreBasicStats[b,2] <- sum(zscoreInsMatrix[b,])
        # Total signal per bp
        zscoreBasicStats[b,3] <- zscoreBasicStats[b,2] / (500 + motifWidth)
        # Motif signal per bp
        zscoreBasicStats[b,4] <- sum(zscoreInsMatrix[b,(250:(250+motifWidth))] / motifWidth)
        # Flank signal per bp
        zscoreBasicStats[b,5] <- (sum(zscoreInsMatrix[b,200:250]) + sum(zscoreInsMatrix[b,(250+motifWidth):(300+motifWidth)])) / 100
        # Background signal per bp
        zscoreBasicStats[b,6] <- (sum(zscoreInsMatrix[b,1:50]) + sum(zscoreInsMatrix[b,(500+motifWidth-50):(500+motifWidth)])) / 100
        # Wide flank signal per bp
        zscoreBasicStats[b,7] <- (sum(zscoreInsMatrix[b,1:250]) + sum(zscoreInsMatrix[b,(250+motifWidth):(500+motifWidth)])) / 500
        # Flank / background
        zscoreBasicStats[b,8] <- zscoreBasicStats[b,5] / zscoreBasicStats[b,6]
        # Motif / flank
        zscoreBasicStats[b,9] <- zscoreBasicStats[b,4] / zscoreBasicStats[b,5]
        # Motif / wide flank
        zscoreBasicStats[b,10] <- zscoreBasicStats[b,4] / zscoreBasicStats[b,7]
      } # end for (b in 1:numSites)
      
      #### Generate null models, use BF and BH correction to parse ####
      
      ## Find the unique values for total signal and generate null models
      uniqueTotalSignals <- unique(siteBasicStats[,2])
      siteWidth <- 500 + motifWidth
      
      ## Initiate a matrix to store the mean null signal in the null model and the input signal to null model
      nullModels <- matrix(data = NA, ncol = 2, nrow = length(uniqueTotalSignals))
      colnames(nullModels) <- c("Input signal", "Avg motif signal in null model")
      
      ## Calculate the null models
      for (c in 1:length(uniqueTotalSignals)){
        nullVec <- generateNullFP(1000, uniqueTotalSignals[c], siteWidth, motifWidth)
        nullModels[c,1] <- uniqueTotalSignals[c]
        nullModels[c,2] <- mean(nullVec)
      } # end for (c in 1:length(uniqueTotalSignals))
      
      ## Perform a one-tailed t-test to generate a p-value for each observed motif site
      cat("Performing one-tailed t-tests on peak subset...", "\n")
      ttest <- list() # list to store the results of the t-tests
      pvalue <- c() # vector to store the p-values
      tvalue <- c() # vector to store the t-value
      
      ## Perform t-test on all sites
      cat("Performing 1-tailed ttest", "\n")
      for (d in 1:numSites){
        ## Retrieve the total signal for the current site
        currentSignal <- c(siteBasicStats[d,2])
        ## Retrieve the appropriate null model
        currentNullModel <- nullModels[which(nullModels[,1]==currentSignal),2]
        ## Perform the t-test
        ttest[[d]] <- t.test(insMatrix[d,250:(250+motifWidth)], mu=currentNullModel, alternative="less", conf.level = 0.95)
        pvalue[d] <- ttest[[d]][["p.value"]]
        tvalue[d] <- ttest[[d]][["statistic"]][["t"]]
      } # end for (d in 1:numSites)
      
      ## Get the indices of the sites that are lower than p = 0.05
      cat("Selecting p-value passing sites > 0.05", "\n")
      idxPvaluePass <- which(pvalue < 0.05)
      pvaluePass <- pvalue[idxPvaluePass]
      ppassNumSites <- length(idxPvaluePass)
      
      ## Perform bonferroni correction
      cat("Performing bonferroni correction", "\n")
      idxBFpass <- which(pvalue < (0.05/numSites))
      BFpvaluePass <- pvalue[idxBFpass]
      BFpassNumSites <- length(idxBFpass)
      
      ## Perform benjamini-hochberg correction
      BHpvalue <- p.adjust(pvalue, method = "BH")
      
      ## Data transfer to storage object and save
      parseData <- list()
      ##
      parseData$PWM <- PWM
      parseData$motifWidth <- motifWidth
      parseData$allSites <- allSites
      parseData$extSites <- extSites
      parseData$numSites <- numSites
      parseData$insMatrix <- insMatrix
      parseData$libSize <- libSize
      parseData$coverageSize <- coverageSize
      parseData$libFactor <- libFactor
      parseData$sampleTotalReads <- sampleTotalReads
      parseData$siteBasicStats <- siteBasicStats
      parseData$samplePeaks <- samplePeaks
      parseData$insStandardDeviation <- insStandardDeviation
      parseData$insMean <- insMean
      parseData$zscoreInsMatrix <- zscoreInsMatrix
      ##
      parseData$uniqueTotalSignals <- uniqueTotalSignals
      parseData$nullModels <- nullModels
      parseData$ttest <- ttest
      parseData$pvalue <- pvalue
      parseData$tvalue <- tvalue
      ##
      parseData$idxPvaluePass <- idxPvaluePass
      parseData$pvaluePass <- pvaluePass
      parseData$ppassNumSites <- ppassNumSites
      ##
      parseData$idxBFpass <- idxBFpass
      parseData$BFpvaluePass <- BFpvaluePass
      parseData$BFpassNumSites <- BFpassNumSites
      ##
      parseData$BHpvalue <- BHpvalue
      ##
      com <- paste0("footprintData$motif", a, "<- parseData")
      eval(parse(text = com))
      
    } # end for (a in 1:numMotif)
    
    ## Save the parsed footprint data
    cat("Saving finished data for", geneName, "\n")
    save(footprintData, file = dataOutPath)
    
  } # end if (numMotif == 0)
} # end if (file.exists(dataOutPath) == TRUE)

## Create the output file for snakemake
file.create(outPath)
cat("Finished parsing", "\n")



