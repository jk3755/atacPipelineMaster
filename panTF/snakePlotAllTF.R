
## Setup : Load Packages #######################################################################################################
library(GenomicRanges)
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
  cat("Processing input file ", a, "\n")
  
  ## Load the footprintData object
  load(fileList[a])
  ## Create a temporary object for the current data
  tempData <- footprintData
  ## Number of potential motifs to plot
  numMotifs <- length(tempData)
  
  if (numMotifs == 1){
    mergedSites <- footprintData[["motif1"]][["peakSites"]]
    mergedRawFootprintMetrics <- footprintData[["motif1"]][["rawFootprintMetrics"]]
  } else {
    
    ## Pull the Granges and insertion matrices from the footprintData object
    for (z in 1:numMotifs){
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
    for (b in 2:numMotifs){
      tryCatch({
        com <- paste0("overlaps <- findOverlaps(sites1, sites", b,")")
        eval(parse(text = com))
        if (length(overlaps@from) == 0){
          # If no overlaps are present, can just directly merge the two Granges
          com <- paste0("sites1 <- c(sites1, sites", b, ")")
          eval(parse(text = com))
          com <- paste0("rawFootprintMetrics1 <- rbind(rawFootprintMetrics1, rawFootprintMetrics", b, ")")
          eval(parse(text = com))
        } else {
          # If overlaps are present, merge only the non-overlapping ranges from the second Granges
          mergeIdx <- overlaps@to
          com <- paste0("sites1 <- c(sites1, sites", b, "[-mergeIdx])")
          eval(parse(text = com))
          com <- paste0("rawFootprintMetrics1 <- rbind(rawFootprintMetrics1, rawFootprintMetrics", b, "[-mergeIdx,])")
          eval(parse(text = com))
        } # end if (length(overlaps@from) == 0)
      }, # end try
      error=function(cond){
        message(cond)
        return(NA)},
      finally={})
    } # end for (b in 2:numMotifs)
    
    mergedSites <- sites1
    mergedRawFootprintMetrics <- rawFootprintMetrics1
    
  } # end if (numMotifs == 1)
  
  ## Try to grab the data for each potential motif
  tryCatch({
  
  ## Calculate the 10% trimmed mean of all insertions in the motif sites
  motifSignal <- mean(mergedRawFootprintMetrics[,3], trim = 0.10)
  ## Calculate the mean of all insertions in the flank region
  flankSignal <- mean(mergedRawFootprintMetrics[,2])
  ## Calculate the mean of background insertions
  backgroundSignal <- mean(mergedRawFootprintMetrics[,1])
  
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
  
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
} # end for (a in 1:numFiles)

### Put the data into dataframe format for plotting

## Count the total number of data points to plot (Total unique motifs)
numPoints <- length(tempBackground)

## get the gene names
x <- substring(fileList, 42)
x <- substring(x, 2)
x <- gsub(".parsedFootprintData.Rdata", "", x)

## Add gene names to dataframe also

## Transfer data to data frame
dfFootprints <- data.frame(flank = tempFlank, depth = tempMotif, background = tempBackground, geneName = x)

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

#################### Code testing


snu61_exp <- coadCounts[,24]

### We need to map between Entrez IDs in the regulon and gene symbols in the other data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene", version = "3.8")
library(mygene)

##

geneList <- 
  
  
queryMany(geneList, scopes="symbol", fields="entrezgene", species="human")
