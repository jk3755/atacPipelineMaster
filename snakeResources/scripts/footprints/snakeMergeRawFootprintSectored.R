
##
cat("Loading libraries", "\n")
suppressMessages(library(GenomicRanges))

##
cat("Setting snakemake vars", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
#
outPath <- snakemake@output[[1]]
#
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Make merged path and perform filecheck
mergedDataPath <- paste0(dirPath, "/footprints/data/merged/", sampleName, "-REP", sampleRep, ".", geneName, ".merged.rawFootprintData.Rdata")
cat("Checking for file at:", mergedDataPath, "\n")

##
if (file.exists(mergedDataPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  ## Scan for all available input files
  sectorPath <- paste0(dirPath, "/footprints/data/temp")
  cat("Scanning for sector footprint files at path:", sectorPath, "\n")
  sectorFiles <- Sys.glob(file.path(sectorPath, paste0(sampleName, "-REP", sampleRep, ".", geneName, "*rawFootprintData.Rdata")))
  cat("Found sector footprint files:", sectorFiles, "\n")
  
  ## Load the files and store in a list
  numFiles <- length(sectorFiles)
  cat("Loading", numFiles, "sectored files", "\n")
  
  ## 
  footprintDataList <- list()
  ##
  for (a in 1:numFiles){
    ##
    load(sectorFiles[a])
    ##
    footprintDataList[[a]] <- footprintData
  } # end for (a in 1:numFiles)
  
  ## Determine the number of motifs
  numMotifs <- length(footprintDataList[[1]])
  cat("Found", numMotifs, "motifs", "\n")
  
  #### Merge the sectored data
  cat("Merging sectored data", "\n")
  
  ##
  footprintData <- list()
  
  ##
  for (b in 1:numMotifs){
    
    ####
    tempList <- list()
    
    #### PWM
    com <- paste0("tempList$PWM <- footprintDataList[[1]][['motif", b, "']][['PWM']]")
    eval(parse(text = com))
    
    #### Motif Width
    com <- paste0("tempList$motifWidth <- footprintDataList[[1]][['motif", b, "']][['motifWidth']]")
    eval(parse(text = com))
    
    #### allSites
    com <- paste0("mergedSites <- c(footprintDataList[[1]][['motif", b, "']][['allSites']], footprintDataList[[2]][['motif", b, "']][['allSites']],
                footprintDataList[[3]][['motif", b, "']][['allSites']], footprintDataList[[4]][['motif", b, "']][['allSites']],
                footprintDataList[[5]][['motif", b, "']][['allSites']], footprintDataList[[6]][['motif", b, "']][['allSites']],
                footprintDataList[[7]][['motif", b, "']][['allSites']], footprintDataList[[8]][['motif", b, "']][['allSites']],
                footprintDataList[[9]][['motif", b, "']][['allSites']], footprintDataList[[10]][['motif", b, "']][['allSites']],
                footprintDataList[[11]][['motif", b, "']][['allSites']], footprintDataList[[12]][['motif", b, "']][['allSites']],
                footprintDataList[[13]][['motif", b, "']][['allSites']], footprintDataList[[14]][['motif", b, "']][['allSites']],
                footprintDataList[[15]][['motif", b, "']][['allSites']], footprintDataList[[16]][['motif", b, "']][['allSites']],
                footprintDataList[[17]][['motif", b, "']][['allSites']], footprintDataList[[18]][['motif", b, "']][['allSites']],
                footprintDataList[[19]][['motif", b, "']][['allSites']], footprintDataList[[20]][['motif", b, "']][['allSites']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$allSites <- mergedSites")
    eval(parse(text = com))
    
    #### numSites
    com <- paste0("tempList$numSites <- length(mergedSites@ranges)")
    eval(parse(text = com))
    
    #### extSites
    com <- paste0("mergedExtSites <- c(footprintDataList[[1]][['motif", b, "']][['extSites']], footprintDataList[[2]][['motif", b, "']][['extSites']],
                footprintDataList[[3]][['motif", b, "']][['extSites']], footprintDataList[[4]][['motif", b, "']][['extSites']],
                footprintDataList[[5]][['motif", b, "']][['extSites']], footprintDataList[[6]][['motif", b, "']][['extSites']],
                footprintDataList[[7]][['motif", b, "']][['extSites']], footprintDataList[[8]][['motif", b, "']][['extSites']],
                footprintDataList[[9]][['motif", b, "']][['extSites']], footprintDataList[[10]][['motif", b, "']][['extSites']],
                footprintDataList[[11]][['motif", b, "']][['extSites']], footprintDataList[[12]][['motif", b, "']][['extSites']],
                footprintDataList[[13]][['motif", b, "']][['extSites']], footprintDataList[[14]][['motif", b, "']][['extSites']],
                footprintDataList[[15]][['motif", b, "']][['extSites']], footprintDataList[[16]][['motif", b, "']][['extSites']],
                footprintDataList[[17]][['motif", b, "']][['extSites']], footprintDataList[[18]][['motif", b, "']][['extSites']],
                footprintDataList[[19]][['motif", b, "']][['extSites']], footprintDataList[[20]][['motif", b, "']][['extSites']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$extSites <- mergedExtSites")
    eval(parse(text = com))
    
    #### insMatrix
    com <- paste0("tempInsMatrix <- rbind(footprintDataList[[1]][['motif", b, "']][['insMatrix']], footprintDataList[[2]][['motif", b, "']][['insMatrix']],
                footprintDataList[[3]][['motif", b, "']][['insMatrix']], footprintDataList[[4]][['motif", b, "']][['insMatrix']],
                footprintDataList[[5]][['motif", b, "']][['insMatrix']], footprintDataList[[6]][['motif", b, "']][['insMatrix']],
                footprintDataList[[7]][['motif", b, "']][['insMatrix']], footprintDataList[[8]][['motif", b, "']][['insMatrix']],
                footprintDataList[[9]][['motif", b, "']][['insMatrix']], footprintDataList[[10]][['motif", b, "']][['insMatrix']],
                footprintDataList[[11]][['motif", b, "']][['insMatrix']], footprintDataList[[12]][['motif", b, "']][['insMatrix']],
                footprintDataList[[13]][['motif", b, "']][['insMatrix']], footprintDataList[[14]][['motif", b, "']][['insMatrix']],
                footprintDataList[[15]][['motif", b, "']][['insMatrix']], footprintDataList[[16]][['motif", b, "']][['insMatrix']],
                footprintDataList[[17]][['motif", b, "']][['insMatrix']], footprintDataList[[18]][['motif", b, "']][['insMatrix']],
                footprintDataList[[19]][['motif", b, "']][['insMatrix']], footprintDataList[[20]][['motif", b, "']][['insMatrix']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$insMatrix <- tempInsMatrix")
    eval(parse(text = com))
    
    #### libSize
    com <- paste0("tempLibSize <- c(footprintDataList[[1]][['motif", b, "']][['libSize']], footprintDataList[[2]][['motif", b, "']][['libSize']],
                footprintDataList[[3]][['motif", b, "']][['libSize']], footprintDataList[[4]][['motif", b, "']][['libSize']],
                footprintDataList[[5]][['motif", b, "']][['libSize']], footprintDataList[[6]][['motif", b, "']][['libSize']],
                footprintDataList[[7]][['motif", b, "']][['libSize']], footprintDataList[[8]][['motif", b, "']][['libSize']],
                footprintDataList[[9]][['motif", b, "']][['libSize']], footprintDataList[[10]][['motif", b, "']][['libSize']],
                footprintDataList[[11]][['motif", b, "']][['libSize']], footprintDataList[[12]][['motif", b, "']][['libSize']],
                footprintDataList[[13]][['motif", b, "']][['libSize']], footprintDataList[[14]][['motif", b, "']][['libSize']],
                footprintDataList[[15]][['motif", b, "']][['libSize']], footprintDataList[[16]][['motif", b, "']][['libSize']],
                footprintDataList[[17]][['motif", b, "']][['libSize']], footprintDataList[[18]][['motif", b, "']][['libSize']],
                footprintDataList[[19]][['motif", b, "']][['libSize']], footprintDataList[[20]][['motif", b, "']][['libSize']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$libSize <- tempLibSize")
    eval(parse(text = com))
    
    #### Coverage Size
    com <- paste0("tempCoverageSize <- c(footprintDataList[[1]][['motif", b, "']][['coverageSize']], footprintDataList[[2]][['motif", b, "']][['coverageSize']],
                footprintDataList[[3]][['motif", b, "']][['coverageSize']], footprintDataList[[4]][['motif", b, "']][['coverageSize']],
                footprintDataList[[5]][['motif", b, "']][['coverageSize']], footprintDataList[[6]][['motif", b, "']][['coverageSize']],
                footprintDataList[[7]][['motif", b, "']][['coverageSize']], footprintDataList[[8]][['motif", b, "']][['coverageSize']],
                footprintDataList[[9]][['motif", b, "']][['coverageSize']], footprintDataList[[10]][['motif", b, "']][['coverageSize']],
                footprintDataList[[11]][['motif", b, "']][['coverageSize']], footprintDataList[[12]][['motif", b, "']][['coverageSize']],
                footprintDataList[[13]][['motif", b, "']][['coverageSize']], footprintDataList[[14]][['motif", b, "']][['coverageSize']],
                footprintDataList[[15]][['motif", b, "']][['coverageSize']], footprintDataList[[16]][['motif", b, "']][['coverageSize']],
                footprintDataList[[17]][['motif", b, "']][['coverageSize']], footprintDataList[[18]][['motif", b, "']][['coverageSize']],
                footprintDataList[[19]][['motif", b, "']][['coverageSize']], footprintDataList[[20]][['motif", b, "']][['coverageSize']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$coverageSize <- tempCoverageSize")
    eval(parse(text = com))
    
    #### libFactor
    com <- paste0("tempLibFactor <- c(footprintDataList[[1]][['motif", b, "']][['libFactor']], footprintDataList[[2]][['motif", b, "']][['libFactor']],
                footprintDataList[[3]][['motif", b, "']][['libFactor']], footprintDataList[[4]][['motif", b, "']][['libFactor']],
                footprintDataList[[5]][['motif", b, "']][['libFactor']], footprintDataList[[6]][['motif", b, "']][['libFactor']],
                footprintDataList[[7]][['motif", b, "']][['libFactor']], footprintDataList[[8]][['motif", b, "']][['libFactor']],
                footprintDataList[[9]][['motif", b, "']][['libFactor']], footprintDataList[[10]][['motif", b, "']][['libFactor']],
                footprintDataList[[11]][['motif", b, "']][['libFactor']], footprintDataList[[12]][['motif", b, "']][['libFactor']],
                footprintDataList[[13]][['motif", b, "']][['libFactor']], footprintDataList[[14]][['motif", b, "']][['libFactor']],
                footprintDataList[[15]][['motif", b, "']][['libFactor']], footprintDataList[[16]][['motif", b, "']][['libFactor']],
                footprintDataList[[17]][['motif", b, "']][['libFactor']], footprintDataList[[18]][['motif", b, "']][['libFactor']],
                footprintDataList[[19]][['motif", b, "']][['libFactor']], footprintDataList[[20]][['motif", b, "']][['libFactor']])")
    eval(parse(text = com))
    ##
    com <- paste0("tempList$libFactor <- tempLibFactor")
    eval(parse(text = com))
    
    #### Transfer merged data to footprintData object
    com <- paste0("footprintData$motif", b, " <- tempList")
    eval(parse(text = com))
    
  } # end for (b in 1:numMotifs)
} # end if (file.exists(mergedDataPath) == TRUE)

#### Output the merged footprintData object
cat("Finished merging, saving data", "\n")
save(footprintData, file = mergedDataPath)


cat("Finished merging!", "\n")
file.create(outPath)








