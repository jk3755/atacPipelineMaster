

## Load libraries
suppressMessages(library(GenomicRanges))

## Set snakemake variables
inputDir <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
dirPath <- snakemake@wildcards[["path"]]

#### Load and process the data from all genes ####
## List all the files that will be analyzed
fileList <- list.files(inputDir, full.names = TRUE)
numFiles <- length(fileList)

cat("Found", numFiles, "footprint data files. Processing...", "\n")

## Initialize a data frame that will store the aggregated data for all TFs
## Will make plot generation much easier
aggregateFootprintData <- data.frame(matrix(vector(), 0, 23,
                          dimnames=list(c(), c(
                          "Gene", "numMotifs", "numPeakSites", "numBoundSites", "numUnboundSites", "boundRatio",
                          "peakMotifSignal", "peakFlankSignal", "peakBackgroundSignal", "peak.log2Flank", "peak.log2Depth",
                          "boundMotifSignal", "boundFlankSignal", "boundBackgroundSignal", "bound.log2Flank", "bound.log2Depth",
                          "unboundMotifSignal", "unboundFlankSignal", "unboundBackgroundSignal", "unbound.log2Flank", "unbound.log2Depth",
                          "RNAexp", "viperNES"))),
                          stringsAsFactors=F)


## Also initialize a list to store the Footprint metrics and site data for each gene
aggregateFootprintMetricsData <- list()

## Use an interation index in case any individual step fails
dataIdx <- 1


## Iterate over all files, aggregate data
for (a in 1:numFiles){
  
  ##
  tryCatch({
  
  load(fileList[a])
    
  ## Use regex to get the gene name of the current file
  geneName <- gsub(inputDir, "", fileList[a])
  geneName <- gsub(paste0(sampleName, "."), "", geneName)
  geneName <- gsub("/", "", geneName)
  geneName <- gsub(".processedFootprintData.Rdata", "", geneName)
  
  #### TRANSFER DATA TO AGGREGATOR OBJECTs ####
  aggregateFootprintData[dataIdx,1] <- geneName
  aggregateFootprintData[dataIdx,2] <- processedFootprintData[["numMotifs"]]
  aggregateFootprintData[dataIdx,3] <- processedFootprintData[["numPeakSites"]]
  aggregateFootprintData[dataIdx,4] <- processedFootprintData[["numBoundSites"]]
  aggregateFootprintData[dataIdx,5] <- processedFootprintData[["numUnboundSites"]]
  aggregateFootprintData[dataIdx,6] <- processedFootprintData[["peakMotifSignal"]]
  aggregateFootprintData[dataIdx,7] <- processedFootprintData[["peakFlankSignal"]]
  aggregateFootprintData[dataIdx,8] <- processedFootprintData[["peakBackgroundSignal"]]
  aggregateFootprintData[dataIdx,9] <- processedFootprintData[["peak.log2Flank"]]
  aggregateFootprintData[dataIdx,10] <- processedFootprintData[["peak.log2Depth"]]
  aggregateFootprintData[dataIdx,11] <- processedFootprintData[["boundMotifSignal"]]
  aggregateFootprintData[dataIdx,12] <- processedFootprintData[["boundFlankSignal"]]
  aggregateFootprintData[dataIdx,13] <- processedFootprintData[["boundBackgroundSignal"]]
  aggregateFootprintData[dataIdx,14] <- processedFootprintData[["bound.log2Flank"]]
  aggregateFootprintData[dataIdx,15] <- processedFootprintData[["bound.log2Depth"]]
  aggregateFootprintData[dataIdx,16] <- processedFootprintData[["unboundMotifSignal"]]
  aggregateFootprintData[dataIdx,17] <- processedFootprintData[["unboundFlankSignal"]]
  aggregateFootprintData[dataIdx,18] <- processedFootprintData[["unboundBackgroundSignal"]]
  aggregateFootprintData[dataIdx,19] <- processedFootprintData[["unbound.log2Flank"]]
  aggregateFootprintData[dataIdx,20] <- processedFootprintData[["unbound.log2Depth"]]
  
  ##
  tempList <- list()
  tempList$"peakSites" <- processedFootprintData[["peakSites"]]
  tempList$"rawPeakFootprintMetrics" <- processedFootprintData[["rawPeakFootprintMetrics"]]
  tempList$"boundSites" <- processedFootprintData[["boundSites"]]
  tempList$"boundSitesMetrics" <- processedFootprintData[["boundSitesMetrics"]]
  tempList$"unboundSites" <- processedFootprintData[["unboundSites"]]
  tempList$"unboundSitesMetrics" <- processedFootprintData[["unboundSitesMetrics"]]
  
  ##
  com <- paste0("aggregateFootprintMetricsData$'", geneName, "' <- tempList")
  eval(parse(text = com))
  
  ##
  dataIdx <- (dataIdx + 1)
  
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
  
} # end for (a in 1:numFiles)

## Combine to one list to save in a single Rdata object
aggregateData <- list()
aggregateData$"aggregateFootprintData" <- aggregateFootprintData
aggregateData$"aggregateFootprintMetricsData" <- aggregateFootprintMetricsData

## Save aggregate data
dataOutPath <- gsub("operations", "data", outPath)
dataOutPath <- gsub("done", "Rdata", dataOutPath)
save(aggregateData, file = dataOutPath)

## Touch the file to update snakemake
file.create(outPath)