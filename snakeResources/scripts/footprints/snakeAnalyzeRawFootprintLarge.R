#### Disable scientific notation in variables
options(scipen = 999)

## Set snakemake variables
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
geneName <- snakemake@wildcards[["gene"]]
currentChunk <- snakemake@wildcards[["chunknum"]]
dirPath <- snakemake@wildcards[["path"]]

## Set the output path for Rdata file and perform a filecheck
footprintDataPath <- paste0(dirPath, "footprints/data/raw/temp/", sampleName, "-REP", sampleRep, ".", geneName, ".chunk", currentChunk, ".rawFootprintData.Rdata")
cat("Output path for raw footprint data:", footprintDataPath, "\n")

##
if (file.exists(footprintDataPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  ## Load libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(stats4))
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(parallel))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(genomation))
  
  ## Load the binding sites data
  cat("Loading binding sites", "\n")
  load(sitesPath)
  numMotif <- length(bindingSites)
  bamFile <- BamFile(bamPath)
  
  ## Initiate an R object to hold all generated data
  footprintData <- list()
  for (a in 1:numMotif){
    com <- paste0("footprintData$motif", a, " <- list()")
    eval(parse(text = com))
  } # end for (a in 1:numMotif)
  
  ##
  cat("Analyzing footprints for", geneName, "\n")
  cat("Found", numMotif, "unique motifs", "\n")
  
  ## If no motifs are found, skip
  if (numMotif ==0){
    cat("No motifs found. Skipping", "\n")
  } else {
    
    ## Index counter for motif naming
    idxMotif <- 1
    
    ## Begin analysis
    for (b in 1:numMotif){
      
      ## Binding Sites
      cat("Analyzing motif", b, "\n")
      PWM <- bindingSites[[b]][["PWM"]]
      motifWidth <- length(bindingSites[[b]][["PWM"]][1,])
      allSites <- bindingSites[[b]][["sites"]]
      allSites <- keepStandardChromosomes(allSites, pruning.mode="coarse")
      allSites <- trim(allSites)
      numSites <- length(allSites)
      cat("Found", numSites, "motif binding sites", "\n")
      
      ## If no binding sites are found, skip this motif
      if (numSites == 0){
        cat("No motif binding sites found, skipping", "\n")
      } else {
        
        ##
        cat("Setting sites for analysis for current chunk", "\n")
        currentChunk <- as.numeric(currentChunk)
        
        ##
        f <- unique(cut(1:numSites,breaks=20))
        labs <- levels(f)[f]
        lower <- as.numeric( sub("\\((.+),.*", "\\1", labs))
        upper <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs))
        
        ##
        if (currentChunk == 1){
          cat("Current chunk is 1", "\n")
          allSites <- allSites[1:upper[1]]
        } else if (currentChunk > 1 && currentChunk < 20){
          cat("Current chunk is between 2-19", "\n")
          allSites <- allSites[upper[(currentChunk - 1)]:upper[currentChunk]]
        } else if (currentChunk == 20){
          cat("Current chunk is 20", "\n")
          allSites <- allSites[upper[19]:numSites]
        }
        
        ##
        numSites <- length(allSites)
        
        ##
        cat("Processing analysis window for each site", "\n")
        ## extend each range +/- 250 bp from motif edges
        extSites <- promoters(allSites, upstream = 250, downstream = (250 + motifWidth), use.names=TRUE)
        extSites <- keepStandardChromosomes(extSites, pruning.mode="coarse")
        extSites <- trim(extSites)
        
        ## Read in the data from bam file for current ranges
        param <- ScanBamParam(which = extSites)
        
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
        shiftedInsertions <- grMerged
        
        ## Perform the footprint calculations
        ## Convert Tn5 insertions corrected Granges to Rle object
        cat("Generating insertion matrix", "\n")
        insRLE <- coverage(grMerged)
        
        ## Create a views object for the Rle list using the Granges sites data
        insViews <- Views(insRLE, extSites)
        
        ## Convert to a matrix
        insMatrix <- as.matrix(insViews)
        
        ## Store the data
        cat("Storing data", "\n")
        ## Initiate a temporary list object to store data, will be transferred to footprintData list
        tempData <- list()
        tempData$PWM <- bindingSites[[b]][["PWM"]]
        tempData$motifWidth <- motifWidth
        tempData$allSites <- allSites
        tempData$numSites <- numSites
        tempData$extSites <- extSites
        tempData$insMatrix <- insMatrix
        tempData$libSize <- length(bamIn)
        tempData$coverageSize <- sum(as.numeric(width(reduce(grIn, ignore.strand=TRUE))))
        tempData$libFactor <- tempData$libSize / tempData$coverageSize
        
        #### Transfer all the data for the current motif to the storage object
        cat("Transferring all data to storage object footprintData", "\n")
        com <- paste0("footprintData$motif", idxMotif, " <- tempData")
        eval(parse(text = com))
        
        ## Update the motif index
        idxMotif <- (idxMotif + 1)
        
      } # end if (numSites == 0)
    } # end for (b in 1:numMotif)
    
    ## Save the raw footprint data
    cat("Saving finished data for", geneName, "\n")
    save(footprintData, file = footprintDataPath)
    
  } # end if (numMotif ==0)
} # end if (file.exists(footprintDataPath) == TRUE)

## Create the output file for snakemake
file.create(outPath)