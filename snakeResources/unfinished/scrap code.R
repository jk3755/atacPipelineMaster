
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
#######################################################################################################

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

