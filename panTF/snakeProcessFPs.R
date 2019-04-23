
## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)

## Disable scientific notation in variables and turn off warnings globally
options(scipen = 999)
options(warn = -1)

## Load libraries
cat("Loading libraries...", "\n")
suppressMessages(library(GenomicRanges))

## Set snakemake variables
cat("Setting snakemake variables...", "\n")
fpDirPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]
##

################################################################################################################################

#### Load and process the data from all genes ####
fileList <- list.files(fpDirPath, full.names = TRUE)
numFiles <- length(fileList)
cat("Found", numFiles, "footprint data files. Processing...", "\n")

#### Iterate over all identified files in the directory and perform the processing for each one ####
for (a in 1:numFiles){
  cat("Processing input file ", a, "\n")
  
  ## Load the current footprintData object
  load(fileList[a])

  #### MERGE AND DEDUPLICATE ALL BINDING SITES ####
  ## At this point, the footprints have already been parsed into bound and unbound sites
  ## Transfer data for both to the new storage object, for downstream analysis
  numMotifs <- length(footprintData)
  
  # If there is only one motif available, there is no need to merge and deduplicate the identified sites
  if (numMotifs == 1){
    
    ##
    peakSites <- footprintData[["motif1"]][["peakSites"]]
    peakInsertionMatrix <- footprintData[["motif1"]][["insMatrix"]]
    rawPeakFootprintMetrics <- footprintData[["motif1"]][["rawFootprintMetrics"]]
    
    ## Split the sites into bound and unbound as determined by null model with bonferroni correction
    boundSiteOverlaps <- findOverlaps(footprintData[["motif1"]][["parseData"]][["bfSites"]], footprintData[["motif1"]][["peakSites"]])
    boundSiteIndex <- boundSiteOverlaps@to
    ##
    boundSites <- peakSites[boundSiteIndex]
    boundSitesInsertionMatrix <- peakInsertionMatrix[boundSiteIndex,]
    boundSitesMetrics <- rawPeakFootprintMetrics[boundSiteIndex,]
    numBoundSites <- length(boundSites)
    ##
    unboundSites <- peakSites[-boundSiteIndex]
    unboundSitesInsertionMatrix <- peakInsertionMatrix[-boundSiteIndex,]
    unboundSitesMetrics <- rawPeakFootprintMetrics[-boundSiteIndex,]
    numUnboundSites <- length(unboundSites)
  
  ## If more than one motif is found for the current TF, merge the data and deduplicate any overlapping genomic sites ##
    
  } else {
    
    ## Pull the data for each individual motif and perform the bound/unbound split
    for (z in 1:numMotifs){
      tryCatch({
        
        ## Pull the basic data
        com <- paste0("peakSites", z, " <- footprintData[['motif", z, "']][['peakSites']]")
        eval(parse(text = com))
        com <- paste0("peakInsertionMatrix", z, " <- footprintData[['motif", z, "']][['insMatrix']]")
        eval(parse(text = com))
        com <- paste0("rawPeakFootprintMetrics", z, " <- footprintData[['motif", z, "']][['rawFootprintMetrics']]")
        eval(parse(text = com))
        
        ## Perform the bound/unbound split
        com <- paste0("boundSiteOverlaps", z, " <- findOverlaps(footprintData[['motif", z, "']][['parseData']][['bfSites']], footprintData[['motif", z, "']][['peakSites']])")
        eval(parse(text = com))
        com <- paste0("boundSiteIndex", z, " <- boundSiteOverlaps", z, "@to")
        eval(parse(text = com))
        
        ## Bound sites
        com <- paste0("boundSites", z, " <- peakSites", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("boundSitesInsertionMatrix", z, " <- peakInsertionMatrix", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics", z, " <- rawPeakFootprintMetrics", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("numBoundSites", z, " <- length(boundSites", z, ")")
        eval(parse(text = com))
        
        ## Unbound sites
        com <- paste0("unboundSites", z, " <- peakSites", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("unboundSitesInsertionMatrix", z, " <- peakInsertionMatrix", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics", z, " <- rawPeakFootprintMetrics", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("numUnboundSites", z, " <- length(unboundSites", z, ")")
        eval(parse(text = com))
        
      }, # end try
      error=function(cond){
        message(cond)
        return(NA)},
      finally={})
    } # end for (z in 1:numMotif)
    
    
    #### Perform the merging and deduplication ####
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
  } # end if (numMotifs == 1)
  
  #### Gather data and transfer to processed storage object ####
  ## Initialize a new list object to store the processed data
  processedFootprintData <- list()
  
  ## Transfer relevant data to the new list
  processedFootprintData$"geneName" <- footprintData[["motif1"]][["geneName"]]
  processedFootprintData$"numMotifs" <- numMotifs
  processedFootprintData$"numUniquePeakSites" <- length(sites1)
  
    ## Transfer the sites and raw signal data to the new storage object
    processedFootprintData$"numUniqueSites" <- length(sites1)
    processedFootprintData$"mergedSites" <- sites1
    processedFootprintData$"mergedSites" <- sites1

    mergedRawFootprintMetrics <- rawFootprintMetrics1
    
  
  
  
  
  
    ## Inititate vectors for temporary data storage
    tempFlank <- c()
    tempMotif <- c()
    tempBackground <- c()
    ## Index iterators for vectors
    idxFlank <- 1
    idxMotif <- 1
    idxBackground <- 1
  
  
  
  
  
  
  
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

#### get the gene names TEMP FIX ####
#genex <- substring(fileList, 42)
#genex <- substring(genex, 2)
#genex <- gsub(".parsedFootprintData.Rdata", "", genex)
##
#genex <- substring(fileList, 39)
#genex <- substring(genex, 2)
#genex <- gsub(".rawFootprintData.Rdata", "", genex)
##
genex <- substring(fileList, 40)
genex <- substring(genex, 3)
genex <- gsub(".rawFootprintData.Rdata", "", genex)
## Add gene names to dataframe

load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
coadCounts <- cclecounts[["large_intestine_bat1"]]
#snu61Exp <- coadCounts[,24]
h508Exp <- coadCounts[,7]
mdst8Exp <- coadCounts[,29]
##
curExp <- mdst8Exp
##
geneList <- names(curExp)
mappings <- queryMany(geneList, scopes="entrezgene", fields="symbol", species="human")
##
genemaps <- mappings@listData[["symbol"]]
##
idx <- which(genemaps %in% genex)


##
addme <- c()
for (m in 1:numPoints){
  f <- which(genemaps %in% genex[m])
  if (length(f)==0){
    addme[m] <- 0
  } else {
    addme[m] <- curExp[[f]]
  }}


##
#for (k in 1:236){
#  idx <- which(n %in% item)
#  val <- curExp[[idx]]
#  dfFootprints$exp[[idx]] <- val}

## Transfer data to data frame
dfFootprints <- data.frame(flank = tempFlank, depth = tempMotif, background = tempBackground, geneName = genex, exp = addme)

#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

####################
load("~/git/atacPipelineMaster/atacVIPER/regulons/coad-tcga-regulon.rda")
#load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
#coadCounts <- cclecounts[["large_intestine_bat1"]]
##
x <- msviper(curExp, regul)
##
plot(x, cex=.7)
##
nex <- x[["es"]][["nes"]]
##
nexname <- names(nex)
nexmappings <- queryMany(nexname, scopes="entrezgene", fields="symbol", species="human")
nexgenemaps <- nexmappings@listData[["symbol"]]
##
idx <- which(nexgenemaps %in% genex)


## Collect the viper values
addmex <- c()
for (m in 1:numPoints){
  f <- which(nexgenemaps %in% genex[m])
  if (length(f)==0){
    addmex[m] <- -10
  } else {
    addmex[m] <- nex[[f]]
  }}


## Transfer data to data frame
dfFootprints$nex <- addmex


## plot the graph with viper values
ggplot(dfFootprints, aes(depth, flank, color=nex)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")



## plot footprint depth against viper activity
ggplot(dfFootprints, aes(nex, depth, color=flank)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")



#########################################################################################################################
# Do the same for bf parsed sites

#### Load and process the data from all genes ####
## List all the files that will be analyzed
fileList <- list.files("C:\\Users\\jsk33\\Desktop\\h508p\\", full.names = TRUE)
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
    footprintData[["motif1"]][["parseData"]][["bfSites"]]
    mergedRawFootprintMetrics <- footprintData[["motif1"]][["parseData"]][["bfFootprintMetrics"]]
  } else {
    
    ## Pull the Granges and insertion matrices from the footprintData object
    for (z in 1:numMotifs){
      tryCatch({
        com <- paste0("sites", z, " <- footprintData[['motif", z, "']][['parseData']][['bfSites']]")
        eval(parse(text = com))
        com <- paste0("rawFootprintMetrics", z, " <- footprintData[['motif", z, "']][['parseData']][['bfFootprintMetrics']]")
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

#### get the gene names TEMP FIX ####
#genex <- substring(fileList, 42)
#genex <- substring(genex, 2)
#genex <- gsub(".parsedFootprintData.Rdata", "", genex)
genex <- substring(fileList, 41)
genex <- substring(genex, 2)
genex <- gsub(".parsedFootprintData.Rdata", "", genex)


## Add gene names to dataframe

load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
coadCounts <- cclecounts[["large_intestine_bat1"]]
#snu61Exp <- coadCounts[,24]
h508Exp <- coadCounts[,7]
##
curExp <- h508Exp
##

geneList <- names(curExp)
mappings <- queryMany(geneList, scopes="entrezgene", fields="symbol", species="human")
##
genemaps <- mappings@listData[["symbol"]]
##
idx <- which(genemaps %in% genex)
##
addme <- c()
for (m in 1:a){
  f <- which(genemaps %in% genex[m])
  if (length(f)==0){
    addme[m] <- 0
  } else {
    addme[m] <- curExp[[f]]
  }}
##
#for (k in 1:nums){
#  idx <- which(n %in% item)
#  val <- snu61exp[[idx]]
#  dfFootprints$exp[[idx]] <- val}

## Transfer data to data frame
dfFootprints <- data.frame(flank = tempFlank, depth = tempMotif, background = tempBackground, geneName = genex, exp = addme )

#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

####################
#### Load and process the data from all genes ####
## List all the files that will be analyzed
#fileList <- list.files("C:\\Users\\jsk33\\Desktop\\h508\\", full.names = TRUE)
fileList <- list.files("C:\\Users\\jsk33\\Desktop\\mdst8\\", full.names = TRUE)
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

#### get the gene names TEMP FIX ####
#genex <- substring(fileList, 42)
#genex <- substring(genex, 2)
#genex <- gsub(".parsedFootprintData.Rdata", "", genex)
##
#genex <- substring(fileList, 39)
#genex <- substring(genex, 2)
#genex <- gsub(".rawFootprintData.Rdata", "", genex)
##
genex <- substring(fileList, 40)
genex <- substring(genex, 3)
genex <- gsub(".rawFootprintData.Rdata", "", genex)
## Add gene names to dataframe

load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
coadCounts <- cclecounts[["large_intestine_bat1"]]
#snu61Exp <- coadCounts[,24]
h508Exp <- coadCounts[,7]
mdst8Exp <- coadCounts[,29]
##
curExp <- mdst8Exp
##
geneList <- names(curExp)
mappings <- queryMany(geneList, scopes="entrezgene", fields="symbol", species="human")
##
genemaps <- mappings@listData[["symbol"]]
##
idx <- which(genemaps %in% genex)


##
addme <- c()
for (m in 1:numPoints){
  f <- which(genemaps %in% genex[m])
  if (length(f)==0){
    addme[m] <- 0
  } else {
    addme[m] <- curExp[[f]]
  }}


##
#for (k in 1:236){
#  idx <- which(n %in% item)
#  val <- curExp[[idx]]
#  dfFootprints$exp[[idx]] <- val}

## Transfer data to data frame
dfFootprints <- data.frame(flank = tempFlank, depth = tempMotif, background = tempBackground, geneName = genex, exp = addme)

#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

####################
load("~/git/atacPipelineMaster/atacVIPER/regulons/coad-tcga-regulon.rda")
#load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
#coadCounts <- cclecounts[["large_intestine_bat1"]]
##
x <- msviper(curExp, regul)
##
#plot(x, cex=.7)
##
nex <- x[["es"]][["nes"]]
##
nexname <- names(nex)
nexmappings <- queryMany(nexname, scopes="entrezgene", fields="symbol", species="human")
nexgenemaps <- nexmappings@listData[["symbol"]]
##
idx <- which(nexgenemaps %in% genex)


## Collect the viper values
addmex <- c()
for (m in 1:numPoints){
  f <- which(nexgenemaps %in% genex[m])
  if (length(f)==0){
    addmex[m] <- -10
  } else {
    addmex[m] <- nex[[f]]
  }}


## Transfer data to data frame
dfFootprints$nex <- addmex


## plot the graph with viper values
ggplot(dfFootprints, aes(depth, flank, color=nex)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")



## plot footprint depth against viper activity
ggplot(dfFootprints, aes(nex, depth, color=flank)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")