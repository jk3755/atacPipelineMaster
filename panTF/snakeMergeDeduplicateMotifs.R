##
library(GenomicRanges)

## Get the number of unique motifs for the current gene
numMotif <- length(footprintData)

## Pull the Granges and insertion matrices from the footprintData object
for (z in 1:numMotif){
  tryCatch({
    com <- paste0("sites", z, " <- footprintData[['motif", z, "']][['peakSites']]")
    eval(parse(text = com))
    com <- paste0("rawFootprintMetrics", z, " <- footprintData[['motif", z, "']][['rawFootprintMetrics']]")
    eval(parse(text = com))
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)},
  finally={})
} # end for (z in 1:numMotif)

## Merge the Granges objects
for (a in 2:numMotif){
  tryCatch({
  com <- paste0("overlaps <- findOverlaps(sites1, sites", a,")")
  eval(parse(text = com))
  if (length(overlaps@from) == 0){
    # If no overlaps are present, can just directly merge the two Granges
    com <- paste0("sites1 <- c(sites1, sites", a, ")")
    eval(parse(text = com))
    com <- paste0("rawFootprintMetrics1 <- rbind(rawFootprintMetrics1, rawFootprintMetrics", a, ")")
    eval(parse(text = com))
  } else {
    # If overlaps are present, merge only the non-overlapping ranges from the second Granges
    mergeIdx <- x@to
    com <- paste0("sites1 <- c(sites1, sites", a, "[-mergeIdx])")
    eval(parse(text = com))
    com <- paste0("rawFootprintMetrics1 <- rbind(rawFootprintMetrics1, rawFootprintMetrics", a, "[-mergeIdx,])")
    eval(parse(text = com))
  } # end if (length(overlaps@from) == 0)
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)},
  finally={})
} # end for (a in 1:numMotif)
##
mergedSites <- sites1
mergedRawFootprintMetrics <- rawFootprintMetrics1

