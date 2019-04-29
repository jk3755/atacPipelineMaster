## Set snakemake variables
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
peakPath <- snakemake@input[[4]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]
##
cat("Processing raw footprint data for gene", geneName, "from sample", sampleName, "\n")

## Set the output path for Rdata file and perform a filecheck
footprintDataPath <- paste0(dirPath, "footprints/data/genome/raw/", sampleName, ".", geneName, ".rawFootprintData.Rdata")
cat("Output path for data is set as:", footprintDataPath, "\n")

## Perform a filecheck first, skip if already exists
if (file.exists(footprintDataPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  ## The file doesn't exist yet, begin analysis
  ## Load libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(genomation))
  suppressMessages(library(stats4))
  suppressMessages(library(rlist))
  
  ## Load the binding sites for current gene
  cat("Loading binding sites", "\n")
  load(sitesPath)
  numMotif <- length(bindingSites)
  bamFile <- BamFile(bamPath)
  
  ## Initiate an R object to hold all generated data
  ## and set a motif number index
  footprintData <- list()
  idxMotif <- 1
  
  ## Loop through all the unique motifs and perform the analysis
  cat("Analyzing footprints for", geneName, "with", numMotif, "unique motifs",  "\n")
  for (b in 1:numMotif){
    
    ## Initiate a temporary list object to store data, will be transferred to footprintData list
    tempData <- list()
    cat("Analyzing motif", b, "\n")
    
    ## Binding Sites
    cat("Processing binding sites", "\n")
    scope <- paste0("chr", c(1:22, "X", "Y"))
    genomeSites <- bindingSites[[b]][["sites"]]
    ## Trim the matched binding sites to the standard chromosomes only
    genomeSites <- keepStandardChromosomes(genomeSites, pruning.mode="coarse")
    genomeSites <- keepSeqlevels(genomeSites, scope, pruning.mode="coarse")
    genomeSites <- trim(genomeSites, use.names = TRUE)
    numSites <- length(genomeSites)
    cat("Found", numSites, "genome-wide binding sites", "\n")
    
    ## Pull binding motif data
    cat("Loading binding motif data", "\n")
    PWM <- bindingSites[[b]][["PWM"]]
    motifWidth <- length(bindingSites[[b]][["PWM"]][1,])
    
    ## Extend each binding site in the Granges to the analysis window (+/- 250 bp)
    cat("Processing analysis window for each binding site", "\n")
    extendedSites <- promoters(genomeSites, upstream = 250, downstream = (250 + motifWidth), use.names=TRUE)
    
    ## Use GenomicAlignments package to read in bam file to GRanges, also very fast
    ## Consider each read as a unique element (insertion), not paired end
    cat("Loading relevant reads", "\n")
    param <- ScanBamParam(which = extendedSites)
    bamIn <- readGAlignments(bamFile, param = param)
      
    ## Convert GAlignments to GRanges
    cat("Converting reads to insertions", "\n")
    grIn <- granges(bamIn)
      
    ## Trim everything but standard chromosomes, trim out of bounds ranges
    cat("Trimming out of bounds reads and keeping only standard chromosomes", "\n")
    grIn <- keepStandardChromosomes(grIn, pruning.mode="coarse")
    grIn <- keepSeqlevels(grIn, scope, pruning.mode="coarse")
    grIn <- trim(grIn, use.names = TRUE)
      
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
      ## Get rid of the mitochondrial data
      insRLE@listData <- insRLE@listData[which(names(insRLE@listData) != "chrM")]
      
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
      
      ## Transfer the data
      tempData$PWM <- bindingSites[[b]][["PWM"]]
      tempData$genomeSites <- allSites
      tempData$motifWidth <- length(bindingSites[[b]][["PWM"]][1,])
      
      tempData$extSites <- promoters(tempData$genomeSites, upstream = 250, downstream = (250 + tempData$motifWidth), use.names=TRUE)
      
      
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
      
      ##
      tempData$rawFootprintMetrics <- rawFootprintMetrics
      
      
      for (a in 1:numMotif){
        com <- paste0("footprintData$motif", a, " <- list()")
        eval(parse(text = com))} # end for (a in 1:numMotif)
      
      #### Transfer all the data for the current motif to the storage object
      cat("Transferring all data to storage object footprintData", "\n")
      com <- paste0("footprintData$motif", idxMotif, " <- tempData")
      eval(parse(text = com))
      
      idxMotif <- (idxMotif + 1)
      
  } # end for (b in 1:numMotif)
  
  ## To avoid errors, clear the list of any empty sub-lists first
  ## Should this result in an object with no data, that can be output
  ## as a dummy file to keep the pipeline running smoothly
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  
  ## Save the raw footprint data
  cat("Saving finished data for", geneName, "\n")
  save(footprintData, file = footprintDataPath)
  
} # end if (file.exists(footprintDataPath) == TRUE)

# Display warnings to the terminal
warnings()

##
file.create(outPath)