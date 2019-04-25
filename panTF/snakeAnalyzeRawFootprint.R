
## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("stats4", suppressUpdates = TRUE)
#biocLite("BiocGenerics", suppressUpdates = TRUE)
#biocLite("parallel", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)
#biocLite("GenomicAlignments", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)

## Disable scientific notation in variables
options(scipen = 999)

## Load libraries
cat("Loading libraries...", "\n")
suppressMessages(library(GenomicRanges))
suppressMessages(library(stats4))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(genomation))

## Set snakemake variables
cat("Setting snakemake variables...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
peakPath <- snakemake@input[[4]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Set the output path for Rdata file and perform a filecheck
footprintDataPath <- paste0(dirPath, "footprints/data/raw/", sampleName, ".", geneName, ".rawFootprintData.Rdata")

if (file.exists(footprintDataPath) == TRUE){
  
  cat("File already exists, skipping", "\n")
  
} else {
  
  ##
  cat("Loading binding sites...", "\n")
  load(sitesPath)
  numMotif <- length(bindingSites)
  bamFile <- BamFile(bamPath)
  
  ## Peaks
  cat("Loading accessibility peaks...", "\n")
  grPeaks <- readBed(peakPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
  grPeaks <- keepStandardChromosomes(grPeaks, pruning.mode="coarse")
  
  ## Initiate an R object to hold all generated data
  footprintData <- list()
  for (a in 1:numMotif){
    com <- paste0("footprintData$motif", a, " <- list()")
    eval(parse(text = com))} # end for (a in 1:numMotif)
  
  cat("Analyzing footprints for", geneName, "\n")
  cat("Found", numMotif, "unique motifs", "\n")
  
  ## Index counter for motif naming, required in case some motifs have no matches in peak sites
  idxMotif <- 1
  
  ## Begin analysis
  for (b in 1:numMotif){
    
    ##
    cat("Analyzing motif", b, "\n")
    
    ## Initiate a temporary list object to store data, will be transferred to footprintData list
    tempData <- list()
    
    cat("Processing binding sites", "\n")
    
    ## Binding Sites
    cat("Subsetting binding sites based on accessibility peaks", "\n")
    allSites <- bindingSites[[b]][["sites"]]
    peakSites <- subsetByOverlaps(allSites, grPeaks)
    numPeakSites <- length(peakSites)
    cat("Found", numPeakSites, "motif binding sites in peak accessibility regions", "\n")
    
    if (numPeakSites == 0){
      
      next
      
    } else {
      
      ## Transfer the data
      tempData$PWM <- bindingSites[[b]][["PWM"]]
      tempData$peakSites <- peakSites
      tempData$numPeakSites <- numPeakSites
      tempData$motifWidth <- length(bindingSites[[b]][["PWM"]][1,])
      
      cat("Processing analysis window for each site", "\n")
      ## extend each range +/- 250 bp from motif edges
      tempData$extSites <- promoters(tempData$peakSites, upstream = 250, downstream = (250 + tempData$motifWidth), use.names=TRUE)
      ## Read in the data from bam file for current ranges
      param <- ScanBamParam(which = tempData$extSites)
      ## Use GenomicAlignments package to read in bam file to GRanges, also very fast
      ## Consider each read as a unique element (insertion), not paired end
      cat("Loading relevant reads", "\n")
      bamIn <- readGAlignments(bamFile, param = param)
      ## Convert GAlignments to GRanges
      cat("Converting reads to insertions", "\n")
      grIn <- granges(bamIn)
      ## Trim everything but standard chromosomes, trim out of bounds ranges
      grIn <- keepStandardChromosomes(grIn, pruning.mode="coarse")
      grIn <- trim(grIn)
      ## Convert the reads to insertions with width = 1
      grIn2 <- resize(grIn, width = 1)
      ## Subset the Granges object into plus and minus strands for shifting
      cat("Shifting insertions +4/-5 bp", "\n")
      plusIdx <- which(strand(grIn2) == "+")
      minusIdx <- which(strand(grIn2) == "-")
      grPlus <- grIn2[plusIdx]
      grMinus <- grIn2[minusIdx]
      ## Shift the ATACseq reads to account for Tn5 insertion mechanism 
      ## shift end of fragment +4 bp (plus strand) or -5 bp (minus standed)
      grPlusShifted <- shift(grPlus, shift=4L)
      grMinusShifted <- shift(grMinus, shift=-5L)
      ## Merge the plus and minus strand shifted Granges
      grMerged <- c(grPlusShifted, grMinusShifted)
      tempData$shiftedInsertions <- grMerged
      
      ## Perform the footprint calculations
      ## Convert Tn5 insertions corrected Granges to Rle object
      cat("Generating insertion matrix", "\n")
      insRLE <- coverage(grMerged)
      ## Get the matching sites
      extSites <- tempData$extSites
      extSites <- keepStandardChromosomes(extSites, pruning.mode="coarse")
      extSites <- trim(extSites)
      ## Create a views object for the Rle list using the Granges sites data
      insViews <- Views(insRLE, extSites)
      ## Convert to a matrix
      insMatrix <- as.matrix(insViews)
      
      ## Calculate the insertion probability at each basepair
      cat("Calculating insertion probabilies", "\n")
      rawTotalSignal <- sum(insMatrix)
      rawProfile <- matrix(data = NA, ncol = length(insMatrix[1,]), nrow = 2)
      rownames(rawProfile) <- c("Column sums", "Insertion frequency")
      ##
      for (c in 1:length(insMatrix[1,])){
        rawProfile[1,c] <- sum(insMatrix[,c])
        rawProfile[2,c] <- (rawProfile[1,c] / rawTotalSignal) * 100
      } # end for (c in 1:length(insMatrix[1,]))
      
      ## Store the data
      cat("Storing data", "\n")
      tempData$extSites <- extSites
      tempData$insRLE <- insRLE
      tempData$insViews <- insViews
      tempData$insMatrix <- insMatrix
      tempData$rawTotalSignal <- rawTotalSignal
      tempData$rawProfile <- rawProfile
      tempData$libSize <- length(bamIn)
      tempData$coverageSize <- sum(as.numeric(width(reduce(grIn, ignore.strand=TRUE))))
      tempData$libFactor <- tempData$libSize / tempData$coverageSize
      ##
      rm(extSites, insRLE, insViews, insMatrix, rawTotalSignal, rawProfile, bamIn)
      rm(grIn, grIn2, plusIdx, minusIdx, grPlus, grMinus, grPlusShifted, grMinusShifted)
      gc()
      
      ## Calculate flanking accessibility and footprint depth data
      cat("Calculating flanking accessibility and footprint depth data", "\n")
      rawFootprintMetrics <- matrix(data = NA, ncol = 5, nrow = length(tempData$insMatrix[,1]))
      colnames(rawFootprintMetrics) <- c("Background", "Flanking", "Motif", "Flanking Accessibility", "Footprint Depth")
      
      for (d in 1:length(tempData$insMatrix[,1])){
        rawFootprintMetrics[d,1] <- (sum(tempData$insMatrix[d,1:50]) + sum(tempData$insMatrix[d,(450 + tempData$motifWidth):(500 + tempData$motifWidth)]))
        rawFootprintMetrics[d,2] <- (sum(tempData$insMatrix[d,200:250]) + sum(tempData$insMatrix[d,(200 + tempData$motifWidth):(250 + tempData$motifWidth)]))
        rawFootprintMetrics[d,3] <- sum(tempData$insMatrix[d,(250:(250 + tempData$motifWidth))])
        rawFootprintMetrics[d,4] <- rawFootprintMetrics[d,2] / rawFootprintMetrics[d,1]
        rawFootprintMetrics[d,5] <- rawFootprintMetrics[d,3] / rawFootprintMetrics[d,2]
      } # end (for d in 1:length(tempData$insMatrix[,1]))
      
      tempData$rawFootprintMetrics <- rawFootprintMetrics
      rm(rawFootprintMetrics)
      gc()
      
      #### Transfer all the data for the current motif to the storage object
      cat("Transferring all data to storage object footprintData", "\n")
      com <- paste0("footprintData$motif", idxMotif, " <- tempData")
      eval(parse(text = com))
      
      ## Update the motif index
      idxMotif <- (idxMotif + 1)
      
    } # end if (numPeakSites = 0)
    
  } # end for (b in 1:numMotif)
  
  ## Save the raw footprint data
  cat("Saving finished data for", geneName, "\n")
  save(footprintData, file = footprintDataPath)
  
} # end if (file.exists(footprintDataPath) == TRUE)

##
file.create(outPath)