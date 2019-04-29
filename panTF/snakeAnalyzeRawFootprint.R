## Set snakemake variables
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sitesPath <- snakemake@input[[3]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]
##
cat("Processing raw footprint data for gene", geneName, "from sample", sampleName, "\n")

#### REMOVE ME #####
sitesPath <- "C:/Users/jsk33/Desktop/bug/MAFF.bindingSites.Rdata"
geneName <- "TEST"
sampleName <- "STEST"
bamPath <- "C:\\Users\\jsk33\\Desktop\\bug\\H508A-WT-02-repmerged.bam"
baiPath <- "C:\\Users\\jsk33\\Desktop\\bug\\H508A-WT-02-repmerged.bai"
#### REMOVE ME #####


## Set the output path for Rdata file and perform a filecheck
footprintDataPath <- paste0(dirPath, "footprints/data/raw/", sampleName, ".", geneName, ".rawFootprintData.Rdata")
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
  
  ## Remove motifs with 0 binding sites
  cat("Removing motifs that matched 0 genomic loci", "\n")
  numMotif <- length(bindingSites)
  numSites <- c()
  ##
  for (l in 1:numMotif){numSites[l] <- length(bindingSites[[l]][["sites"]]@ranges)}
  zeroIdx <- which(numSites == 0)
  bindingSites <- bindingSites[-zeroIdx]
  numMotif <- length(bindingSites)
  
  ## Subset to unique PWM only ##
  cat("Removing duplicate motifs", "\n")
  uniqueBindingSites <- list()
  uniqueBindingSites[1] <- bindingSites[1]
  PWMidx <- 2
  ##
  for (x in 2:numMotif){
    addPWM <- "YES"
    curNumPWM <- length(uniqueBindingSites)
    curBindingSites <- bindingSites[[x]][["sites"]]
    ##
    for (z in 1:curNumPWM){
      compBindingSites <- uniqueBindingSites[[z]][["sites"]]
      cat(addPWM, "\n")
      ##
      if (identical(curBindingSites, compBindingSites)){
        cat("sites are identical", "\n")
        addPWM <- "NO"}} # end for (z in 1:curNumPWM)}
    
    if (addPWM == "YES"){
      uniqueBindingSites[PWMidx] <- bindingSites[x]
      PWMidx <- (PWMidx + 1)}
  } # end for (x in 2:numMotif)
  ##
  bindingSites <- uniqueBindingSites
  
  ## Initiate an R object to hold all generated data
  ## and set a motif number index
  bamFile <- BamFile(bamPath)
  footprintData <- list()
  idxMotif <- 1
  
  ## Loop through all the unique motifs and perform the analysis
  cat("Analyzing footprints for", geneName, "with", numMotif, "unique motifs",  "\n")
  for (b in 1:numMotif){
    
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
      
    ## Using a tryCatch block, errors in any given motif won't stop pipeline
    #tryCatch({
    
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
    ## Convert reads to insertions
    grIn <- resize(grIn, width = 1)

    ## Subset the Granges object into plus and minus strands for shifting
    ## Shift the ATACseq reads to account for Tn5 insertion mechanism 
    ## shift end of fragment +4 bp (plus strand) or -5 bp (minus standed)
    cat("Shifting insertions by +4/-5 bp", "\n")
    plusIdx <- which(strand(grIn) == "+")
    minusIdx <- which(strand(grIn) == "-")
    grPlus <- grIn[plusIdx]
    grMinus <- grIn[minusIdx]
    grPlusShifted <- shift(grPlus, shift=4L)
    grMinusShifted <- shift(grMinus, shift=-5L)
    ## Merge the plus and minus strand shifted Granges
    grShiftedInsertions <- c(grPlusShifted, grMinusShifted)
    
    ## Perform the footprint calculations
    ## Convert Tn5 insertions corrected Granges to Rle object
    cat("Generating insertion matrix", "\n")
    insertionRLE <- coverage(grShiftedInsertions)
    ## Get rid of the mitochondrial data
    insertionRLE@listData <- insertionRLE@listData[which(names(insertionRLE@listData) != "chrM")]
    ## Get the matching sites
    extendedSites <- keepStandardChromosomes(extendedSites, pruning.mode="coarse")
    extendedSites <- keepSeqlevels(extendedSites, scope, pruning.mode="coarse")
    extendedSites <- trim(extendedSites, use.names = TRUE)
    ## Create a views object for the Rle list using the Granges sites data
    insertionViews <- Views(insertionRLE, extendedSites)
    ## Convert to a matrix
    insertionMatrix <- as.matrix(insertionViews)
      
    ## Store and save all the data for downstream analysis
    cat("Storing data in a list object", "\n")
    tempData <- list()
    ## Transfer the data
    tempData$librarySize <- length(bamIn)
    tempData$coverageSize <- sum(as.numeric(width(reduce(grIn, ignore.strand=TRUE))))
    tempData$libraryFactor <- (tempData$librarySize / tempData$coverageSize)
    ##
    tempData$PWM <- PWM
    tempData$genomeSites <- genomeSites
    tempData$numGenomeSites <- numSites
    tempData$motifWidth <- motifWidth
    tempData$extendedSites <- extendedSites
    tempData$shiftedInsertions <- grShiftedInsertions
    tempData$insertionRLE <- insertionRLE
    tempData$insertionViews <- insertionViews
    tempData$insertionMatrix <- insertionMatrix

    ## Transfer all the data for the current motif to the storage object
    cat("Transferring all data to storage object footprintData", "\n")
    com <- paste0("footprintData$motif", idxMotif, " <- tempData")
    eval(parse(text = com))
    
    ## Update the motif index counter
    cat("Updating motif index counter", "\n")
    idxMotif <- (idxMotif + 1)
    
    ## Cleanup variables
    rm(tempData, insertionMatrix, insertionViews, extendedSites, insertionRLE, grShiftedInsertions,
       grMinusShifted, grPlusShifted, grMinus, grPlus, minusIdx, plusIdx, grIn, bamIn, genomeSites)
    gc()
    
    # }, # end try
    # error = function(cond){
    #   message(cond)
    #   return(NA)
    # },
    # finally={})
    
   
    
  } # end for (b in 1:numMotif)
  
  ## Save the raw footprint data
  cat("Saving finished raw footprint data for", geneName, "\n")
  save(footprintData, file = footprintDataPath)
  
} # end if (file.exists(footprintDataPath) == TRUE)

## Create the output file for snakemake
cat("Creating output file for snakemake", geneName, "\n")
file.create(outPath)
gc()
warnings()