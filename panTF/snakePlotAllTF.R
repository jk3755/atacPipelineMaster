
## Setup : Load Packages #######################################################################################################
loadLibrary <- function(lib){
  # Require returns TRUE invisibly if it was able to load package
  # If package was not able to be loaded then first check if it can be loaded with Bioconductor
  # If not, attempt to reinstall with Bioconductor and then install.packages 
  # Load package after (re)installing
  for(i in lib){
    if (!require("BiocManager")){
      install.packages("BiocManager")
        if(!require(i, character.only = TRUE)){
          BiocManager::install(i)
            if(!require(i, character.only = TRUE)){
              install.packages(i, dependencies = TRUE)
                require(i, character.only = TRUE)}}}}}
## Required libraries
packages <- c(
              "ggplot2",
              "ggsci")

loadLibrary(packages)
options(warn = -1)
options(scipen = 999)
################################################################################################################################


## Plot 1 - Footprint depth against flanking accessibility, for all TFs, dots colored by RNA expression level

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
  tempData <- footprintData
  
  ## Number of potential motifs to plot
  numMotifs <- length(tempData)
  
  ## Try to grab the data for each potential motif
  tryCatch({
  
  ##
  com <- paste0("tempMotifData <- tempData$motif", a)
  eval(parse(text = com))
  ##
  motifWidth <- tempMotifData[["motifWidth"]]
  tempMetrics <- tempMotifData[["rawFootprintMetrics"]]
  
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
  idxFlank <- 1
  idxMotif <- 1
  idxBackground <- 1


  
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
  
  
  tempMotifData <- tempData[1]
  
  ##
  dataMatrix[a,1] <- mean(temp[,4])
  dataMatrix[a,2] <- mean(temp[,5])
  dataMatrix[a,3] <- mean(temp[,1])
  
}


## Count the total number of data points to plot (Number of total unique motifs for all unique TFs)
##
dataMatrix <- matrix(data = NA, nrow = numFiles, ncol = 3)
colnames(dataMatrix) <- c("Flank", "Motif", "Background")
##
data <- data.frame(flank = dataMatrix[,1], footprint = dataMatrix[,2], background = dataMatrix[,3])

##
ggplot(data, aes(footprint, flank, color=background)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

##
ggplot(data, aes(footprint, flank, color=background)) + 
  geom_point() + 
  scale_fill_gsea()



## Plot 2 - Footprint depth against flanking accessibility, for all TFs, dots colored by VIPER protein activity level
