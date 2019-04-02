
## Setup : Load Packages #######################################################################################################
library(ggplot2)
library(ggsci)
##
options(warn = -1)
options(scipen = 999)
################################################################################################################################

#### Load and process the data from all genes ####
## List all the files that will be analyzed
fileList <- list.files("C:\\Users\\jsk33\\Desktop\\fpdata\\", full.names = TRUE)
numFiles <- length(fileList)

## Inititate vectors for temporary data storage
tempFlank <- c()
tempMotif <- c()
tempBackground <- c()
## Index iterators for vectors
idxFlank <- 1
idxMotif <- 1
idxBackground <- 1

## Iterate over each unique transcription factor
for (a in 1:numFiles){
  
  ## Load the footprintData object
  load(fileList[a])
  ## Create a temporary object for the current data
  tempData <- footprintData
  ## Number of potential motifs to plot
  numMotifs <- length(tempData)
  
  ## Iterate over all potential motifs
  for (b in 1:numMotifs){
  
  ## Try to grab the data for each potential motif
  tryCatch({
  
  com <- paste0("tempMotifData <- tempData$motif", b)
  eval(parse(text = com))
  ##
  motifWidth <- tempMotifData[["motifWidth"]]
  tempMetrics <- tempMotifData[["rawFootprintMetrics"]]
  
  ## Catch an error where the loaded motif list has no data
  if (is.null(tempMetrics)){next}else{
  
  ## Calculate the 10% trimmed mean of all insertions in the motif sites
  motifSignal <- mean(tempMetrics[,3], trim = 0.10)
  ## Calculate the mean of all insertions in the flank region
  flankSignal <- mean(tempMetrics[,2])
  ## Calculate the mean of background insertions
  backgroundSignal <- mean(tempMetrics[,1])
  
  ## Calculate flanking accessibility (log2 fold change between flank and background)
  log2Flank <- log2(flankSignal/backgroundSignal)
  ## Calculate footprint depth (log2 fold change between flank and background)
  log2Depth <- log2(motifSignal/flankSignal)
  
  ## Transfer the data to the temporary vectors
  tempFlank[idxFlank] <- log2Flank
  tempMotif[idxMotif] <- log2Depth
  tempBackground[idxBackground] <- backgroundSignal
  
  ## Update vector indices
  idxFlank <- (idxFlank +1)
  idxMotif <- (idxMotif +1)
  idxBackground <- (idxBackground +1)
  
  } # end if (is.null(tempMetrics))
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
  } # end for (b in 1:numMotifs)
} # end for (a in 1:numFiles)

## Put the data into dataframe format for plotting

## Count the total number of data points to plot (Total unique motifs)
numPoints <- length(tempBackground)
## Transfer data to data frame
dfFootprints <- data.frame(flank = tempFlank, depth = tempMotif, background = tempBackground)

#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=background)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

##
ggplot(data, aes(footprint, flank, color=background)) + 
  geom_point() + 
  scale_fill_gsea()



## Plot 2 - Footprint depth against flanking accessibility, for all TFs, dots colored by VIPER protein activity level
