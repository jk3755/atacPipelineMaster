#### Generate gene names in snakemake group format and output to text file. Only need to run once
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

#### Library loading
cat("Loading libraries...", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))

## Set the output filepath
outPath <- snakemake@output[[1]]

## Query the database to return all human annotations
mdbHuman <- query (MotifDb, 'hsapiens')
## Get the list of all unique genes
uniqueGenes <- unique(mdbHuman@elementMetadata@listData[["geneSymbol"]])
## Get the total number of unique annotated genes
numGenes <- length(uniqueGenes)
## Initiate list object to store all motif data
motifData <- list()

##
for (gene in uniqueGenes){
  
  ## Get relevant indices
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
  ## Number of redundant motifs for current gene
  numMotifs <- length(geneIdx)
  
  ## Get the unique motifs for current gene
  tempMotifs <- list()
  c <- 1
  
  ##
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1} # end for (idx in geneIdx)
  
  ##
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  
  ## Initiate sublists for current gene
  tryCatch({
    com <- paste0("motifData$", gene, " <- list()")
    eval(parse(text = com))},
    error=function(cond){
      #message(cond)
      return(NA)},
    finally={}) #end tryCatch()
  
  ## Populate motifData list for current gene
  tryCatch({     
    for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
  },
  error=function(cond){
    #message(cond)
    return(NA)},
  finally={}) #end tryCatch()
  
} # end for (gene in uniqueGenes)

## Save the motifData object
save(motifData, file = outPath)


#### Generate gene uniqueGenes in snakemake group format and output to text file. Only need to run once
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

#### Library loading
cat("Loading libraries...", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))

## Query the database to return all human annotations
mdbHuman <- query (MotifDb, 'hsapiens')
## Get the list of all unique genes
uniqueGenes <- unique(mdbHuman@elementMetadata@listData[["geneSymbol"]])
## Get the total number of unique annotated genes
numGenes <- length(uniqueGenes)
## Initiate list object to store all motif data
motifData <- list()

##
for (gene in uniqueGenes){
  
  ## Get relevant indices
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
  ## Number of redundant motifs for current gene
  numMotifs <- length(geneIdx)
  
  ## Get the unique motifs for current gene
  tempMotifs <- list()
  c <- 1
  
  ##
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1} # end for (idx in geneIdx)
  
  ##
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  
  ## Initiate sublists for current gene
  tryCatch({
    com <- paste0("motifData$", gene, " <- list()")
    eval(parse(text = com))},
    error=function(cond){
      #message(cond)
      return(NA)},
    finally={}) #end tryCatch()
  
  ## Populate motifData list for current gene
  tryCatch({     
    for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
  },
  error=function(cond){
    #message(cond)
    return(NA)},
  finally={}) #end tryCatch()
  
} # end for (gene in uniqueGenes)

##################### Output text file
## Set the output filepath
outPath <- snakemake@output[[1]]
numGenes <- length(motifData)
uniqueGenes <- names(motifData)
strings <- c()

##
a <- 1 # count for genes
b <- 1 # string index

##
while (a <= numGenes){
  
  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4
  h <- a+5
  i <- a+6
  j <- a+7
  k <- a+8
  l <- a+9
  m <- a+10
  n <- a+11
  o <- a+12
  p <- a+13
  q <- a+14
  r <- a+15
  s <- a+16
  t <- a+17
  u <- a+18
  v <- a+19
  
  tmp1 <- paste0("'sites/operations/", uniqueGenes[c], ".PWMscan.done", "', ")
  tmp2 <- paste0("'sites/operations/", uniqueGenes[d], ".PWMscan.done", "', ")
  tmp3 <- paste0("'sites/operations/", uniqueGenes[e], ".PWMscan.done", "', ")
  tmp4 <- paste0("'sites/operations/", uniqueGenes[f], ".PWMscan.done", "', ")
  tmp5 <- paste0("'sites/operations/", uniqueGenes[g], ".PWMscan.done", "', ")
  tmp6 <- paste0("'sites/operations/", uniqueGenes[h], ".PWMscan.done", "', ")
  tmp7 <- paste0("'sites/operations/", uniqueGenes[i], ".PWMscan.done", "', ")
  tmp8 <- paste0("'sites/operations/", uniqueGenes[j], ".PWMscan.done", "', ")
  tmp9 <- paste0("'sites/operations/", uniqueGenes[k], ".PWMscan.done", "', ")
  tmp10 <- paste0("'sites/operations/", uniqueGenes[l], ".PWMscan.done", "', ")
  tmp11 <- paste0("'sites/operations/", uniqueGenes[m], ".PWMscan.done", "', ")
  tmp12 <- paste0("'sites/operations/", uniqueGenes[n], ".PWMscan.done", "', ")
  tmp13 <- paste0("'sites/operations/", uniqueGenes[o], ".PWMscan.done", "', ")
  tmp14 <- paste0("'sites/operations/", uniqueGenes[p], ".PWMscan.done", "', ")
  tmp15 <- paste0("'sites/operations/", uniqueGenes[q], ".PWMscan.done", "', ")
  tmp16 <- paste0("'sites/operations/", uniqueGenes[r], ".PWMscan.done", "', ")
  tmp17 <- paste0("'sites/operations/", uniqueGenes[s], ".PWMscan.done", "', ")
  tmp18 <- paste0("'sites/operations/", uniqueGenes[t], ".PWMscan.done", "', ")
  tmp19 <- paste0("'sites/operations/", uniqueGenes[u], ".PWMscan.done", "', ")
  tmp20 <- paste0("'sites/operations/", uniqueGenes[v], ".PWMscan.done", "'")
  
  strings[b] <- paste0(
    "rule PWMscan_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'sites/operations/PWMscan.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'\n"
  )
  a <- a+20
  b <- b+1
  
}

## Write the file
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)



## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("stats4", suppressUpdates = TRUE)
#biocLite("BiocGenerics", suppressUpdates = TRUE)
#biocLite("parallel", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)
#biocLite("GenomicAlignments", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)

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
footprintDataPath <- paste0(dirPath, "footprints/data/genome/raw/", sampleName, ".", geneName, ".rawFootprintData.Rdata")

if (file.exists(footprintDataPath) == TRUE){
  
  cat("File already exists, skipping", "\n")
  
} else {
  
  ## Load libraries
  cat("Loading libraries...", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(stats4))
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(parallel))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(genomation))
  suppressMessages(library(rlist))
  
  ##
  cat("Loading binding sites...", "\n")
  load(sitesPath)
  numMotif <- length(bindingSites)
  bamFile <- BamFile(bamPath)
  
  ## Initiate an R object to hold all generated data
  footprintData <- list()
  for (a in 1:numMotif){
    com <- paste0("footprintData$motif", a, " <- list()")
    eval(parse(text = com))} # end for (a in 1:numMotif)
  
  cat("Analyzing footprints for", geneName, "\n")
  cat("Found", numMotif, "unique motifs", "\n")
  
  ##
  idxMotif <- 1
  
  ## Begin analysis
  for (b in 1:numMotif){
    
    ##
    cat("Analyzing motif", b, "\n")
    
    ## Initiate a temporary list object to store data, will be transferred to footprintData list
    tempData <- list()
    
    cat("Processing binding sites", "\n")
    
    ## Binding Sites
    allSites <- bindingSites[[b]][["sites"]]
    ## Trim the matched binding sites to the standard chromosomes only
    scope <- paste0("chr", c(1:22, "X", "Y"))
    allSites <- keepStandardChromosomes(allSites, pruning.mode="coarse")
    allSites <- keepSeqlevels(allSites, scope, pruning.mode="coarse")
    ##
    numSites <- length(allSites)
    cat("Found", numSites, "genome-wide binding sites", "\n")
    
    ## Transfer the data
    tempData$PWM <- bindingSites[[b]][["PWM"]]
    tempData$genomeSites <- allSites
    tempData$motifWidth <- length(bindingSites[[b]][["PWM"]][1,])
    
    cat("Processing analysis window for each site", "\n")
    ## extend each range +/- 250 bp from motif edges
    tempData$extSites <- promoters(tempData$genomeSites, upstream = 250, downstream = (250 + tempData$motifWidth), use.names=TRUE)
    
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
footprintDataPath <- paste0(dirPath, "footprints/data/peaks/raw/", sampleName, ".", geneName, ".rawFootprintData.Rdata")

if (file.exists(footprintDataPath) == TRUE){
  
  cat("File already exists, skipping", "\n")
  
} else {
  
  ## Load libraries
  cat("Loading libraries...", "\n")
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(stats4))
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(parallel))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(genomation))
  suppressMessages(library(rlist))
  
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
    ## Trim everything but standard chromosomes
    scope <- paste0("chr", c(1:22, "X", "Y"))
    allSites <- keepStandardChromosomes(allSites, pruning.mode="coarse")
    allSites <- keepSeqlevels(allSites, scope, pruning.mode="coarse")
    ##
    peakSites <- subsetByOverlaps(allSites, grPeaks)
    numPeakSites <- length(peakSites)
    cat("Found", numPeakSites, "motif binding sites in peak accessibility regions", "\n")
    
    if (numPeakSites == 0){
      
      cat("No binding motifs found in accessible regions for current motif, skipping", "\n")
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


## See below link for reference information related to ChIPseeker package
# https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

## Install libraries, if necessary
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("ChIPseeker")
#BiocManager::install("genomation")
#BiocManager::install("GenomicRanges")
#BiocManager::install("clusterProfiler")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ReactomePA")

##
cat("Loading libraries...", "\n")
library(ChIPseeker)
library(genomation)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

##
cat("Setting snakemake vars...", "\n")
bedFile <- snakemake@input[[1]]
outputPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
sampleRep <- snakemake@wildcards[["repnum"]]

## Data import
cat("Input peak file path:", bedFile, "\n")
narrowPeaks <- ChIPseeker::readPeakFile(bedFile)
## Subset to standard xsomes
narrowPeaks <- keepStandardChromosomes(narrowPeaks, pruning.mode="coarse")


## Coverage plots make plot of genome wide peak coverage
covplotPath <- gsub("annotations.done.txt", "repmerged.peakgenomecov.svg", outputPath)
covplotPath <- gsub("operations", "metrics", covplotPath)
cat("Output path for peak genomve coverage plot:", covplotPath, "\n")
##
cat("Generating genome-wide peak coverage plot...", "\n")
svg(file = covplotPath) # set the filepath for saving the svg figure
##
weightname <- names(narrowPeaks@elementMetadata@listData[2])
covplot(narrowPeaks, weightCol=weightname)
## Turn off svg device 
dev.off()


## Profile of peaks in TSS regions
TSSprofilePath <- gsub("peakgenomecov", "peakTSSprofile", covplotPath)
cat("Output path for peak TSS profile plot:", TSSprofilePath, "\n")
## One step function to generate TSS heatmap from a BED file
cat("Generating peak TSS profile plot...", "\n")
svg(file = TSSprofilePath) # set the filepath for saving the svg figure
##
peakHeatmap(bedFile, TxDb=txdb, upstream=3000, downstream=3000, color="red")
## Turn off svg device 
dev.off()


## Average profile of ChIP peaks binding to TSS region
AvgPeakProfileTSS <- gsub("peakgenomecov", "avgPeakTSSprofile", covplotPath)
cat("Output path for average peak TSS profile plot:", AvgPeakProfileTSS, "\n")
## One step function of above from a BED file
cat("Generating average peak TSS profile plot...", "\n")
svg(file = AvgPeakProfileTSS) # set the filepath for saving the svg figure
##
plotAvgProf2(bedFile, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
## Turn off svg device 
dev.off()


## With confidence interval estimate and bootstrap methods
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)


## Peak annotations
cat("Generating peak annotations...", "\n")
peakAnno <- annotatePeak(bedFile, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
## Plot it
annoPlot1 <- gsub("peakgenomecov", "annoplot1", covplotPath)
annoPlot2 <- gsub("peakgenomecov", "annoplot2", covplotPath)
annoPlot3 <- gsub("peakgenomecov", "annoplot3", covplotPath)
annoPlot4 <- gsub("peakgenomecov", "annoplot4", covplotPath)
#annoPlot5 <- gsub("peakgenomecov", "annoplot1", covplotPath)
##
svg(file = annoPlot1) # set the filepath for saving the svg figure
plotAnnoPie(peakAnno)
dev.off()
##
svg(file = annoPlot2) # set the filepath for saving the svg figure
plotAnnoBar(peakAnno)
dev.off()
##
svg(file = annoPlot3) # set the filepath for saving the svg figure
vennpie(peakAnno)
dev.off()
##
svg(file = annoPlot4) # set the filepath for saving the svg figure
upsetplot(peakAnno)
dev.off()
##
#svg(file = annoPlot5) # set the filepath for saving the svg figure
#upsetplot(peakAnno, vennpie=TRUE)


## Finish up
file.create(outputPath)


##
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#BiocManager::install("viper")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("mygene")
#install.packages("rlist")

##
cat("Loading libraries...", "\n")
suppressMessages(library(aracne.networks))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(annotate))
suppressMessages(library(rlist))
suppressMessages(library(viper))
suppressMessages(library(GenomicRanges))
suppressMessages(library(mygene))

##
cat("Setting snakemake vars...", "\n")
inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
entrezid <- snakemake@wildcards[["entrez"]]

## load file
load(inputfile)

## get binding sites of merged signals
sites <- mergedMotifs[["sites"]]

## get the ARACNe interactome
coad_interactome <- aracne.networks::reguloncoad

## get the targets
com <- paste0("targets <- names(coad_interactome[['", entrezid, "']][['tfmode']])")
eval(parse(text = com))


## Retrieve info for the network targets
target_list <- list()
for (a in 1:length(targets)){
  target_list[a] <- mygene::getGene(geneid = targets[a], fields = "all")}


## count the number of gene locations we have
loc_count <- 0
for (a in 1:length(targets)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      loc_count <- (loc_count + length(target_list[[a]][["genomic_pos"]]))
    } else {
      loc_count <- (loc_count + 1)}}}


## Retrieve the genomic coordinates of all network targets
## retrieve all genomic coords here, can later trim non standard ones more easily with GRanges functions
target_locations <- matrix(data = NA, nrow = loc_count, ncol = 4)
colnames(target_locations) <- c("gene", "chr", "start", "end")
idx <- 1
#
for (a in 1:length(target_list)){
  
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      
      for (b in 1:length(target_list[[a]][["genomic_pos"]])){
        target_locations[idx,1] <- target_list[[a]][["symbol"]]
        target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][[b]][["chr"]])
        target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][[b]][["start"]]
        target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][[b]][["end"]]
        idx <- (idx+1)
      } # end for (b in 1:length(target_list[[a]][["genomic_pos"]]))
      
    } else {
      
      target_locations[idx,1] <- target_list[[a]][["symbol"]]
      target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][["chr"]])
      target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][["start"]]
      target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][["end"]]
      idx <- (idx+1)}}}


## Convert the genomic locations of the targets into GRanges
gr <- GRanges(
  seqnames = target_locations[,2],
  ranges = IRanges(start = as.numeric(target_locations[,3]), end = as.numeric(target_locations[,4])))
#strand = c(rep("+", times = num_names)),
#seqinfo = Seqinfo(genome="hg38"),
#score = c(rep(1, times = num_names)))

## prune to standard xsomes
gr <- keepStandardChromosomes(gr, pruning.mode="coarse")

gr1 <- gr
gr2 <- gr
#
start(gr1) <- (start(gr1) - 2000)
width(gr1) <- (width(gr1) + 2000)
#
start(gr2) <- (start(gr2) - 50000)
width(gr2) <- (width(gr2) + 50000)


##
#intersection <- intersect(gr, wg_sites)

## Find overlaps
## YOU WILL WANT TO FIND THE OVERLAPS FOR A RANGE OF EXTENDED VALUES PAST THE TARGETS, GRAPH IT
overlaps1 <- findOverlaps(gr1, sites)
overlaps2 <- findOverlaps(gr2, sites)

info <- list()
info$overlap1 <- overlaps1
info$overlap2 <- overlaps2

save(info, file = outputfile)




cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(MotifDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))

##
cat("Setting snakemake vars...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]
mergedpath <- snakemake@input[[4]]
#
outpathdone <- snakemake@output[[1]]
#
samplename <- snakemake@wildcards[["sample"]]
genename <- snakemake@wildcards[["gene"]]
dirpath <- snakemake@wildcards[["path"]]
prob <- snakemake@wildcards[["prob"]]

##
load(sitespath)
x <- 1 # set to motif 1

graphpath <- paste0(dirpath, "preprocessing/15downsample/footprints/graphs/", samplename, ".", prob, ".", genename, ".svg")
cat("Output path for signal object: ", graphpath, "\n")

if (file.exists(graphpath) == TRUE){
  cat("File already exists, skipping...", "\n")
  next
  
} else {
  
  ##
  cat("Setting parameters...", "\n")
  motif_score <- "99%"
  upstream <- 100
  downstream <- 100
  genome <- Hsapiens
  anchor <- c("cut site")
  
  ##
  cat("Loading data...", "\n")
  mergein <- gsub("done.merged.txt", "merged.Rdata", mergedpath)
  load(mergein)
  sigs <- merged_signals[["signal"]]
  
  ##
  PWM <- bindingSites[[1]][["PWM"]]
  wid <- length(PWM[1,])
  
  cat("Calculating libFactor...", "\n")
  mt <- merged_signals$bindingSites
  scope <- seqlevels(mt)
  cat("Pruning sites 1...", "\n")
  mt <- keepStandardChromosomes(mt, pruning.mode="coarse")
  #
  cat("Pruning sites 2...", "\n")
  mt <- keepSeqlevels(mt, scope, pruning.mode="coarse")
  mt$userdefined <- TRUE
  seqlevels(mt) <- scope
  seqinfo(mt) <- Seqinfo(scope, seqlengths = seqlengths(mt))
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which)
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                    bampath, baipath, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                    bampath, baipath, SIMPLIFY = FALSE)
  }
  ##
  bamIn <- lapply(bamIn, as, Class = "GRanges")
  if(!is(bamIn, "GRangesList")) bamIn <- GRangesList(bamIn)
  bamIn <- unlist(bamIn)
  seqlevelsStyle(bamIn) <- seqlevelsStyle(genome)
  if(anchor=="cut site"){
    ## keep 5'end as cutting sites
    bamIn <- promoters(bamIn, upstream=0, downstream=1)
  }else{
    ## keep fragment center
    bamIn <- reCenterPeaks(bamIn, width=1)
  }
  ##
  libSize <- length(bamIn)
  coverageSize <- sum(as.numeric(width(reduce(bamIn, ignore.strand=TRUE))))
  libFactor <- libSize / coverageSize
  
  ##
  cat("Building profile...", "\n")
  ## segmentation the signals
  ## x2 because stranded.
  Profile <- lapply(sigs, function(.ele) colMeans(.ele, na.rm = TRUE)*2/libFactor)
  ## upstream + wid + downstream
  Profile.split <- lapply(Profile, function(.ele){
    list(upstream=.ele[seq.int(upstream)],
         binding=.ele[upstream+seq.int(wid)],
         downstream=.ele[upstream+wid+seq.int(downstream)])
  })
  
  ##
  optimalSegmentation <- function(.ele){
    .l <- length(.ele)
    short_abun <- cumsum(.ele)/seq.int(.l)
    long_abun <- cumsum(rev(.ele))/seq.int(.l)
    long_abun <- rev(long_abun)
    short_abun <- short_abun[-length(short_abun)]
    long_abun <- long_abun[-1]
    ##long_abun should always greater than short_abun
    long_abun <- long_abun - short_abun
    long_abun[long_abun<0] <- 0
    cov_diff <- numeric(length(short_abun))
    for(i in seq_along(.ele)){
      cov_diff_tmp <- .ele
      cov_diff_tmp <- cov_diff_tmp-short_abun[i]
      cov_diff_tmp[-seq.int(i)] <- cov_diff_tmp[-seq.int(i)] - long_abun[i]
      cov_diff[i] <- mean(cov_diff_tmp^2)
    }
    .ids <- which(cov_diff==min(cov_diff, na.rm = TRUE))
    data.frame(pos=.ids, short_abun=short_abun[.ids], long_abun=long_abun[.ids])
  }
  Profile.seg <- lapply(Profile.split, function(.ele){
    ups <- optimalSegmentation(.ele$upstream)
    downs <- optimalSegmentation(rev(.ele$downstream))
    ## find the nearest pair
    .min <- c(max(rbind(ups, downs)), 0, 0)
    for(i in seq.int(nrow(ups))){
      for(j in seq.int(nrow(downs))){
        tmp <- sum(abs(ups[i, -1] - downs[j, -1]))
        if(tmp < .min[1]){
          .min <- c(tmp, i, j)
        }
      }
    }
    c(colMeans(rbind(ups[.min[2], ], downs[.min[3], ])), binding=mean(.ele$binding, na.rm=TRUE))
  })
  
  ##
  Profile.seg <- colMeans(do.call(rbind, Profile.seg))
  Profile.seg[3] <- Profile.seg[2]+Profile.seg[3]
  names(Profile.seg)[2:3] <- c("distal_abun", "proximal_abun")
  Profile <- c(Profile[["+"]], Profile[["-"]])
  
  ## Make graph
  #svg_path <- "/home/ubuntu1/atac/ls1034/wt01/graphs/test.svg"
  svg(file = graphpath) # set the filepath for saving the svg figure
  cat("Saving svg footprint image at path:", graphpath, "\n")
  
  ## FP graph
  pwm2pfm <- function(pfm, name="motif"){
    if(!all(round(colSums(pfm), digits=4)==1)){
      return(NULL)
    }
    new("pfm", mat=as.matrix(pfm), name=name)
  }
  PWMin <- pwm2pfm(PWM)
  cat("Plotting graph...", "\n")
  ATACseqQC:::plotFootprints(Profile,
                             Mlen=wid, motif=PWMin)
  dev.off()
  
}


##
cat("Finished...", "\n")
file.create(outpathdone)

cat("Loading libraries...", "\n")
#BiocManager::install("ATACseqQC")
suppressMessages(library(ATACseqQC))

##
cat("Setting snakemake vars...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outputSVG <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
#
sampleLabel <- paste0(sampleName, "-", sampleRep)

##
cat("Output filepath for fragment size distribution: ", outputSVG, "\n")
svg(file = outputSVG) # set the filepath for saving the svg figure

## Calculate the fragment size distribution and save to svg
fragSizeDist(
  bamFiles = bamPath,
  bamFiles.labels = sampleLabel,
  index = baiPath
)

## Turn off svg device 
dev.off()

##
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(MotifDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))

##
cat("Setting snakemake vars...", "\n")
final_outpath <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
repnum <- snakemake@wildcards[["repnum"]]
reptot <- snakemake@wildcards[["reptot"]]
genename <- snakemake@wildcards[["gene"]]
current_chr <- snakemake@wildcards[["chr"]]
sample_path <- snakemake@wildcards[["path"]]
current_prob <- snakemake@wildcards[["prob"]]

## Set the input files
cat("Setting input file paths...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]

##
cat("Setting parameters...", "\n")
motif_score <- "99%"
upstream <- 100
downstream <- 100
scope <- current_chr
genome <- Hsapiens

##
cat("Loading binding sites and PWM...", "\n")
load(sitespath)
num_motif <- length(bindingSites)

## Only do motif 1
a <- 1
#
signalpath <- paste0(sample_path, "saturation/footprints/data/bychr/", samplename, "-REP", repnum, "of", reptot, ".", current_prob, ".", genename, ".", current_chr, ".Rdata")
cat("Output path for signal object: ", signalpath, "\n")

#
if (file.exists(signalpath) == TRUE){
  cat("File already exists, skipping...", "\n")
  next
} else {
  
  #
  cat("File not found, proceeding...", "\n")
  PWM <- bindingSites[[a]][["PWM"]]
  sites <- bindingSites[[a]][["sites"]]
  #
  
  cat("Pruning sites 1...", "\n")
  sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
  #
  cat("Pruning sites 2...", "\n")
  sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
  
  ## error handling
  # if current sites object has < 5 sites, the pipeline will crash
  
  #
  cat("Generating signal for ", genename, "motif", a, "chromosome ", current_chr, "\n")
  # generate signal
  sigs <- tryCatch(
    {
      factorFootprints(bamfiles = bampath,
                       index = baipath,
                       bindingSites = sites,
                       pfm = PWM,
                       genome = genome,
                       min.score = motif_score,
                       seqlev = scope,
                       upstream = upstream,
                       downstream = downstream)
      
      # Save the data
      cat("Saving signals data...", "\n")
      save(sigs, file = signalpath)
    },
    error=function(cond){
      message(cond)
      return(NA)
    },
    finally={
    })
} # end if (file.exists(signalpath))

#
cat("Finished...", "\n")
file.create(final_outpath)


## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("stats4", suppressUpdates = TRUE)
#biocLite("BiocGenerics", suppressUpdates = TRUE)
#biocLite("parallel", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)
#biocLite("GenomicAlignments", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
#biocLite("seqLogo", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)
#install.packages("ggplot2")
#install.packages("ggpubr")

## Disable scientific notation in variables
options(scipen = 999)

## Set snakemake variables
cat("Setting snakemake variables", "\n")
footprintDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Set the output filepath for the Rdata object and perform a filecheck
dataOutPath <- gsub("operations", "data", outPath)
dataOutPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", dataOutPath, perl = TRUE)

cat("Output path for parsed data:", dataOutPath, "\n")

if (file.exists(dataOutPath) == TRUE){
  
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
  suppressMessages(library(seqLogo))
  suppressMessages(library(ChIPpeakAnno))
  suppressMessages(library(rlist))
  
  ## Load the footprintData file
  cat("Loading footprintData file", "\n")
  footprintDataPath <- gsub("operations", "data", footprintDataPath)
  footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
  load(footprintDataPath)
  
  ## To avoid errors, clear the list of any empty sub-lists first
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  
  ## If the data object is empty, skip the parse operation and output a dummy file
  if (length(footprintData) == 0){
    
    cat("No data found in footprint object. Skipping", "\n")
    
  } else {
    
    ## Build functions
    cat("Building functions", "\n")
    
    generateNullFP <- function(iterations, inputSignal, analysisWidth, motifWidth){
      # This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
      #- Consider the total signal (number of insertions) at each specific ~200 bp locus
      #- Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
      #- Use published or experimentally derived models of Tn5 sequence specific insertion bias
      #- For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
      #- Generate the null model graph by weighted random residstribution of the total observed signal at that site
      #- Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
      # These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
      # iterations = number of iterations
      # inputSignals = unique values for total signal
      # analysisWidth = total bp in region of interest (flank + background + motif)
      # motifWidth = motif width
      
      ##
      cat("Generating a null footprint model with the following parameters:", "\n")
      cat("Iterations:", iterations, "\n")
      cat("Input signal:", inputSignal, "\n")
      cat("Analysis window (bp):", analysisWidth, "\n")
      cat("Motif width (bp):", motifWidth, "\n")
      
      # declare vector of size n to store average motif signal values
      averages <- c()
      
      # generate the null models and calculate motif averages
      for (a in 1:iterations){
        
        # declare the null vector
        null <- c(1:(analysisWidth))
        
        # randomly distribute the total signal
        # size = the number of values to distribute
        # prob = probability of each site
        # length = length of the generated vector
        null <- c(as.vector(rmultinom(1, size = inputSignal, prob = rep(1, length(null)))))
        
        ## Calculate the mean signal in motif region
        motifStart <- ((analysisWidth - motifWidth) / 2)
        motifEnd <- (motifStart + motifWidth)
        motifAvg <- (sum(null[motifStart:motifEnd])) / motifWidth
        
        ## Store the average values
        averages[a] <- motifAvg
        
      } # end for (a in 1:n)
      return(averages)
    } # end generateNullFP function
    
    plotInsProb <- function(plotTitle = c(""), motifWidth, motifPWM, plotLogo = FALSE, insVector){
      
      ## This function uses code adapted from the R package ATACSeqQC
      ## Plot the figure in a new page in the viewport
      grid.newpage()
      
      ## Data
      totalBP <- length(insVector)
      flankBP <- ((totalBP - motifWidth) / 2 ) ## The number of BP flanking the motif on each side
      
      ## Plot information
      xlab = "Dist. to motif (bp)"
      ylab = "Tn5 fragmentation probability"
      xlim <- c(0, totalBP + 1)
      ylim <- c(0, max(insVector) * 1.12)
      
      ## Add the plotting margins to the viewport (sets the outer bounds of the entire image)
      vp <- plotViewport(margins=c(5.1, 5.1, 4.1, 2.1), name="plotRegion")
      pushViewport(vp)
      
      ## Viewport for the graph plotting area
      vp1 <- viewport(y=.4, height=.8,
                      xscale=xlim,
                      yscale=ylim,
                      name="footprints")
      pushViewport(vp1)
      
      ## Add the insertion probability data line
      grid.lines(x=1:totalBP,
                 y=insVector,
                 default.units="native",
                 gp=gpar(lwd = 1, col = "darkred")) # lwd = line width, col = line color
      
      ## This code adds the x and y axis lines
      # at = is a numeric vector with the x-axis locations for tick marks
      grid.xaxis(at = 
                   c(seq(1, flankBP, length.out = 3),
                     flankBP + seq(1, motifWidth),
                     flankBP + motifWidth + seq(1, flankBP, length.out = 3)),
                 label = c(-(flankBP + 1 - seq(1, flankBP + 1, length.out = 3)),
                           rep("", motifWidth),
                           seq(0, flankBP, len = 3)))
      grid.yaxis()
      
      ## Adds the dashed line across the x-axis horizontally (motif hashes)
      grid.lines(x=c(flankBP, flankBP, 0), y=c(0, max(insVector), ylim[2]),
                 default.units="native", gp=gpar(lty=2))
      
      ##
      grid.lines(x=c(flankBP + motifWidth + 1, flankBP + motifWidth + 1, totalBP),
                 y=c(0, max(insVector), ylim[2]),
                 default.units="native", gp=gpar(lty=2))
      upViewport()
      
      ##
      vp2 <- viewport(y=.9, height=.2,
                      xscale=c(0, totalBP + 1),
                      name="motif")
      pushViewport(vp2)
      upViewport()
      
      ##
      legvp <- viewport(x=0.5,
                        y=0.5,
                        width=convertX(unit(1, "lines"), unitTo="npc"),
                        height=convertY(unit(1, "lines"), unitTo="npc"),
                        just=c("right", "top"), name="legendWraper")
      pushViewport(legvp)
      upViewport()
      
      ##
      grid.text(plotTitle,
                y=unit(1, "npc")-convertY(unit(1, "lines"), unitTo="npc"),
                gp=gpar(cex=1.2, fontface="bold"))
      upViewport()
      
      ## Add the x and y axis labels to the image
      grid.text(xlab, y=unit(1, 'lines'))
      grid.text(ylab, x=unit(1, 'line'), rot = 90)
      
    } # end plotInsProb function
    
    #### PARSING ANALYSIS ####
    
    ## The number of unique motifs for the current gene
    numMotif <- length(footprintData)
    cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
    
    ##
    for (a in 1:numMotif){
      
      ##
      cat("Processing motif", a, "\n")
      
      ## will need to improve this code at some point
      #tryCatch({
      
      ## Prepare the data
      cat("Loading data", "\n")
      com <- paste0("tempData <- footprintData$motif", a)
      eval(parse(text = com))
      genomeSites <- tempData[["genomeSites"]]
      numSites <- length(tempData[["insMatrix"]][,1])
      siteWidth <- length(tempData[["insMatrix"]][1,])
      motifWidth <- tempData[["motifWidth"]]
      PWM <- footprintData[["motif1"]][["PWM"]][a]
      insMatrix <- tempData[["insMatrix"]]
      insProfile <- tempData[["rawProfile"]]
      insVector <- insProfile[2,]
      siteTotalSignal <- c()
      
      ## Make graph of the raw peak sites
      #svgPath <- paste0(dirPath, "footprints/graphs/peaks/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
      #svg(file = svgPath)
      #cat("Saving peaks footprint image at path:", svgPath, "\n")
      #plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".rawpeaks")
      #plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = insVector)
      #dev.off()
      
      ## Calculate total signal for each site
      cat("Calculating total signal at each site", "\n")
      ##
      for (b in 1:numSites){
        siteTotalSignal[b] <- sum(tempData[["insMatrix"]][b,])
      } # end for (b in 1:numSites)
      
      ## Find the unique values for total signal and generate null models
      uniqueTotalSignals <- unique(siteTotalSignal)
      ## Remove NA values from uniqueTotalSignals
      ## (how do they get there???)
      uniqueTotalSignals <- uniqueTotalSignals[!is.na(uniqueTotalSignals)]
      
      ## Initiate a matrix to store the mean null signal in the null model and the input signal to null model
      nullModels <- matrix(data = NA, ncol = 2, nrow = length(uniqueTotalSignals))
      colnames(nullModels) <- c("Input signal", "Avg motif signal in null model")
      
      ## Calculate the null models
      cat("Generating null models", "\n")
      for (c in 1:length(uniqueTotalSignals)){
        nullVec <- generateNullFP(1000, uniqueTotalSignals[c], siteWidth, motifWidth)
        nullModels[c,1] <- uniqueTotalSignals[c]
        nullModels[c,2] <- mean(nullVec)
      } # end for (c in 1:length(uniqueTotalSignals))
      
      ## Perform a one-tailed t-test to generate a p-value for each observed motif site
      cat("Performing one-tailed t-tests on genome-wide binding sites", "\n")
      ttestGenome <- list() # list to store the results of the t-tests
      pvalueGenome <- c() # vector to store the p-values
      tvalueGenome <- c() # vector to store the t-value
      
      ## Perform t-test on all sites
      for (d in 1:numSites){
        ## Retrieve the total signal for the current site
        currentSignal <- c(siteTotalSignal[d])
        ## Retrieve the appropriate null model
        currentNullModel <- nullModels[which(nullModels[,1]==currentSignal),2]
        ## Perform the t-test
        ttestGenome[[d]] <- t.test(insMatrix[d,250:(250+motifWidth)], mu=currentNullModel, alternative="less", conf.level = 0.95)
        pvalueGenome[d] <- ttestGenome[[d]][["p.value"]]
        tvalueGenome[d] <- ttestGenome[[d]][["statistic"]][["t"]]
      } # for (d in 1:numSites)
      
      ## Get the indices of the sites that are lower than p = 0.05
      cat("Selecting p-value passing sites", "\n")
      idxPvaluePass <- which(pvalueGenome < 0.05)
      genomePvaluePass <- pvalueGenome[idxPvaluePass]
      
      ## Perform bonferroni correction
      cat("Performing bonferroni correction", "\n")
      idxbfGenomePass <- which(pvalueGenome < (0.05/numSites))
      bfPvalueGenomePass <- pvalueGenome[idxbfGenomePass]
      
      ## Subset the insertion matrix based on the bonferroni passing sites only
      cat("Subsetting sites based on bf corrected p-values", "\n")
      bfInsMatrix <- insMatrix[idxbfGenomePass,]
      
      ##
      bfTotalSignal <- sum(bfInsMatrix)
      bfProfile <- matrix(data = NA, ncol = length(bfInsMatrix[1,]), nrow = 2)
      rownames(bfProfile) <- c("Column sums", "Insertion frequency")
      
      for (e in 1:length(bfInsMatrix[1,])){
        bfProfile[1,e] <- sum(bfInsMatrix[,e])
        bfProfile[2,e] <- (bfProfile[1,e] / bfTotalSignal) * 100
      } # end for (e in 1:length(bfInsMatrix[1,]))
      
      ##
      bfVector <- bfProfile[2,]
      bfSites <- genomeSites[idxbfGenomePass]
      bfNumSites <- length(idxbfGenomePass)
      
      ## Make a plot for the bf passing sites
      #svgPath <- paste0(dirPath, "footprints/graphs/bf/", sampleName, ".", geneName, ".", "motif", a, ".bf.sites.svg")
      #svg(file = svgPath)
      #cat("Saving peaks footprint image at path:", svgPath, "\n")
      #plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".bfsites")
      #plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = bfVector)
      #dev.off()
      
      # ## HEATMAPS ## CURRENTLY OFF ##
      # #### Make heatmap for bf passing sites ####
      # ## USE THIS STRUCTURE FOR HEATMAPS
      # ## first, combine the signals from plus and minus strand
      # heatSigs <- bfInsMatrix
      # heatNumSites <- bfNumSites
      # heatNumBP <- siteWidth
      # heatSites <- bfSites
      # ## scale each row individually
      # for (f in 1:heatNumSites){
      #   maxsig <- max(heatSigs[f,])
      #   for (g in 1:heatNumBP){heatSigs[f,g] <- (heatSigs[f,g] / maxsig)}}
      # maxsig <- 1
      # ## invert signals
      # for (h in 1:heatNumSites){for (i in 1:heatNumBP){heatSigs[h,i] <- (1-heatSigs[h,i])}}
      # 
      # ## Annotate the combined sublist name which will become the tital of the heatmap plot
      # heatTitle <- paste0(geneName, "_motif", a, "_numsites", heatNumSites)
      # combined <- list()
      # com <- paste0("combined$", heatTitle, " <- heatSigs")
      # eval(parse(text = com))
      # svgPath <- paste0(dirPath, "footprints/graphs/heatmaps/", sampleName, ".", geneName, ".", "motif", a, ".bfpeak.sites.heatmap.svg")
      # svg(file = svgPath)
      # cat("Saving svg footprint image at path:", svgPath, "\n")
      # ## Margin controls
      # # margin(a,b,c,d)
      # # a = size of graph from top to bottom, higher value = smaller. default = 0.1
      # # b = size of graph from left to right, higher value = smaller. default = 0.005
      # # c = flips x axis?
      # # d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
      # # good settings for ATACseq = c(0.1,0.005,0.05,0.2)
      # # bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values
      # ChIPpeakAnno::featureAlignedHeatmap(combined,
      #                                     feature.gr = reCenterPeaks(heatSites, width = heatNumBP),
      #                                     upper.extreme = maxsig, # set this to control the heatmap scale
      #                                     annoMcols = "score",
      #                                     sortBy = "score",
      #                                     n.tile = heatNumBP,
      #                                     margin = c(0.1, 0.005, 0.05, 0.2),
      #                                     color = colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias = 0.9)(100),
      #                                     gp = gpar(fontsize = 10),
      #                                     newpage = TRUE)
      # dev.off()
      
      ## Data transfer to storage object and save
      cat("Transferring data to storage object", "\n")
      parseData <- list()
      ##
      parseData$numSites <- numSites
      parseData$insVector <- insVector
      parseData$siteTotalSignal <- siteTotalSignal
      parseData$uniqueTotalSignals <- uniqueTotalSignals
      parseData$nullModels <- nullModels
      parseData$ttestGenome <- ttestGenome
      parseData$pvalueGenome <- pvalueGenome
      parseData$tvalueGenome <- tvalueGenome
      parseData$genomePvaluePass <- genomePvaluePass
      parseData$bfPvalueGenomePass <- bfPvalueGenomePass
      parseData$bfInsMatrix <- bfInsMatrix
      parseData$bfTotalSignal <- bfTotalSignal
      parseData$bfProfile <- bfProfile
      parseData$bfVector <- bfVector
      parseData$bfSites <- bfSites
      parseData$bfNumSites <- bfNumSites
      ##
      com <- paste0("footprintData$motif", a, "$parseData <- parseData")
      eval(parse(text = com))
      
      #}, # end try
      #error=function(cond){
      #  message(cond)
      #  return(NA)
      #},
      #finally={})
      
    } # end for (a in 1:numMotif)
  } # end if (length(footprintData) == 0)
  
  ## Save data or create dummy file
  save(footprintData, file = dataOutPath)
  
} # end if (file.exists(dataOutPath) == TRUE)

##
gc()
file.create(outPath)
cat("Finished!", "\n")




## Install libraries, if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges", suppressUpdates = TRUE)

## Disable scientific notation in variables and turn off warnings globally
options(scipen = 999)
options(warn = -1)

## Set snakemake variables
cat("Setting snakemake variables", "\n")
inputPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]
bamCopy <- snakemake@wildcards[["bamcopy"]]

## Set the filepaths and perform a filecheck
cat("Setting filepaths and checking for output file", "\n")
inputPath <- gsub("operations", "data", inputPath)
inputPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", inputPath, perl = TRUE)
dataOutPath <- gsub("parsed", "processed", inputPath)

tryCatch({
  
  ##
  if (file.exists(dataOutPath) == TRUE){
    cat("File already exists, skipping", "\n")
  } else {
    
    ## Load libraries
    cat("Loading libraries", "\n")
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(rlist))
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    ## Load the current parsed footprintData object
    cat("Loading parsed footprint file", "\n")
    load(inputPath)
    ## To avoid errors, clear the list of any empty sub-lists first
    cat("Cleaning list", "\n")
    footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
    ## Also remove any sublists that do not have the parsed data object stored
    footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 15L, TRUE)
    numMotifs <- length(footprintData)
    cat("Found", numMotifs, "motifs", "\n")
    
    ## Error handling
    if (numMotifs == 0){
      
      cat("No motifs found. Exiting", "\n")
      
    } else if (numMotifs == 1){
      
      ## Because some motifs may be removed due to errors, pull the motif names to be used in downstream commands (can't just go sequentially)
      motifNames <- names(footprintData)
      
      #### MERGE AND DEDUPLICATE ALL BINDING SITES ####
      ## At this point, the footprints have already been parsed into bound and unbound sites
      ## Transfer data for both to the new storage object, for downstream analysis
      # If there is only one motif available, there is no need to merge and deduplicate the identified sites
      cat("Processing 1 motif", "\n")
      cat("Setting genome sites", "\n")
      com <- paste0("genomeSites <- footprintData[['", motifNames[1], "']][['genomeSites']]")
      eval(parse(text = com))
      cat("Setting genome insertion matrix", "\n")
      com <- paste0("genomeInsertionMatrix <- footprintData[['", motifNames[1], "']][['insMatrix']]")
      eval(parse(text = com))
      cat("Setting footprint metrics", "\n")
      com <- paste0("rawGenomeFootprintMetrics <- footprintData[['", motifNames[1], "']][['rawFootprintMetrics']]")
      eval(parse(text = com))
      numGenomeSites <- length(genomeSites)
      cat("Number of genome sites", numGenomeSites, "\n")
      ## Split the sites into bound and unbound as determined by null model with bonferroni correction
      cat("Finding bound site overlaps", "\n")
      com <- paste0("boundSiteOverlaps <- findOverlaps(footprintData[['", motifNames[1], "']][['parseData']][['bfSites']], footprintData[['", motifNames[1], "']][['genomeSites']])")
      eval(parse(text = com))
      cat("Setting index for bound sites", "\n")
      boundSiteIndex <- boundSiteOverlaps@to
      cat("Subsetting bound sites", "\n")
      boundSites <- genomeSites[boundSiteIndex]
      boundSitesInsertionMatrix <- genomeInsertionMatrix[boundSiteIndex,]
      boundSitesMetrics <- rawGenomeFootprintMetrics[boundSiteIndex,]
      numBoundSites <- length(boundSites)
      cat("Subsetting unbound sites", "\n")
      unboundSites <- genomeSites[-boundSiteIndex]
      unboundSitesInsertionMatrix <- genomeInsertionMatrix[-boundSiteIndex,]
      unboundSitesMetrics <- rawGenomeFootprintMetrics[-boundSiteIndex,]
      numUnboundSites <- length(unboundSites)
      
      #### Calculate footprint characteristics on merged data ####
      ## Calculate the 10% trimmed mean of all insertions in the motif sites
      
      ## I am going to add in code here to convert these values to a per bp average
      ## Should move the code to parse script at some point
      com <- paste0("motifWidth <- footprintData[['", motifNames[1], "']][['motifWidth']]")
      eval(parse(text = com))
      ##
      genomeMotifSignal <- (mean(rawGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
      boundMotifSignal <- (mean(boundSitesMetrics[,3], trim = 0.10) / motifWidth)
      unboundMotifSignal <- (mean(unboundSitesMetrics[,3], trim = 0.10) / motifWidth)
      
      ## Calculate the mean of all insertions in the flank region
      genomeFlankSignal <- (mean(rawGenomeFootprintMetrics[,2], trim = 0.10) / 100)
      boundFlankSignal <- (mean(boundSitesMetrics[,2], trim = 0.10) / 100)
      unboundFlankSignal <- (mean(unboundSitesMetrics[,2], trim = 0.10) / 100)
      
      ## Calculate the mean of background insertions
      genomeBackgroundSignal <- (mean(rawGenomeFootprintMetrics[,1], trim = 0.10) / 100)
      boundBackgroundSignal <- (mean(boundSitesMetrics[,1], trim = 0.10) / 100)
      unboundBackgroundSignal <- (mean(unboundSitesMetrics[,1], trim = 0.10) / 100)
      
      ## Calculate flanking accessibility (log2 fold change between flank and background)
      genome.log2Flank <- log2(genomeFlankSignal / genomeBackgroundSignal)
      bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
      unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
      
      ## Calculate footprint depth (log2 fold change between flank and background)
      genome.log2Depth <- log2(genomeMotifSignal / genomeFlankSignal)
      bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
      unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
      
      #### TRANSFER DATA TO STORAGE OBJECT ####
      ## Initialize a new list object to store the processed data
      processedFootprintData <- list()
      ##
      processedFootprintData$"geneName" <- geneName
      processedFootprintData$"numMotifs" <- numMotifs
      processedFootprintData$"numGenomeSites" <- numGenomeSites
      processedFootprintData$"numBoundSites" <- numBoundSites
      processedFootprintData$"numUnboundSites" <- numUnboundSites
      ##
      processedFootprintData$"genomeSites" <- genomeSites
      processedFootprintData$"rawGenomeFootprintMetrics" <- rawGenomeFootprintMetrics
      processedFootprintData$"genomeMotifSignal" <- genomeMotifSignal
      processedFootprintData$"genomeFlankSignal" <- genomeFlankSignal
      processedFootprintData$"genomeBackgroundSignal" <- genomeBackgroundSignal
      processedFootprintData$"genome.log2Flank" <- genome.log2Flank
      processedFootprintData$"genome.log2Depth" <- genome.log2Depth
      ##
      processedFootprintData$"boundSites" <- boundSites
      processedFootprintData$"boundSitesMetrics" <- boundSitesMetrics
      processedFootprintData$"boundMotifSignal" <- boundMotifSignal
      processedFootprintData$"boundFlankSignal" <- boundFlankSignal
      processedFootprintData$"boundBackgroundSignal" <- boundBackgroundSignal
      processedFootprintData$"bound.log2Flank" <- bound.log2Flank
      processedFootprintData$"bound.log2Depth" <- bound.log2Depth
      ##
      processedFootprintData$"unboundSites" <- unboundSites
      processedFootprintData$"unboundSitesMetrics" <- unboundSitesMetrics
      processedFootprintData$"unboundMotifSignal" <- unboundMotifSignal
      processedFootprintData$"unboundFlankSignal" <- unboundFlankSignal
      processedFootprintData$"unboundBackgroundSignal" <- unboundBackgroundSignal
      processedFootprintData$"unbound.log2Flank" <- unbound.log2Flank
      processedFootprintData$"unbound.log2Depth" <- unbound.log2Depth
      
      
      #### CODE TESTING - PROMOTER/DISTAL groups ####
      ## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
      promoters <- promoters(txdb, upstream = 1000, downstream = 100)
      ## Trim the GRanges object to keep standard entries only
      scope <- paste0("chr", c(1:22, "X", "Y"))
      promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
      promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
      
      ## Subset based on the overlaps
      promoterOverlaps <- findOverlaps(promoters, genomeSites, ignore.strand = TRUE)
      promoterIdx <- unique(promoterOverlaps@to)
      ##
      promoterGenomeSites <- genomeSites[promoterIdx]
      distalGenomeSites <- genomeSites[-promoterIdx]
      ##
      promoterBoundOverlaps <- findOverlaps(promoterGenomeSites, boundSites)
      promoterUnboundOverlaps <- findOverlaps(promoterGenomeSites, unboundSites)
      promoterBoundIdx <- unique(promoterBoundOverlaps@from)
      promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
      ##
      promoterBoundSites <- genomeSites[promoterBoundIdx]
      promoterUnboundSites <- genomeSites[promoterUnboundIdx]
      ##
      distalBoundOverlaps <- findOverlaps(distalGenomeSites, boundSites)
      distalUnboundOverlaps <- findOverlaps(distalGenomeSites, unboundSites)
      distalBoundIdx <- unique(distalBoundOverlaps@from)
      distalUnboundIdx <- unique(distalUnboundOverlaps@from)
      ##
      distalBoundSites <- genomeSites[distalBoundIdx]
      distalUnboundSites <- genomeSites[distalUnboundIdx]
      
      ##### promoter genome sites
      promoterGenomeFootprintMetrics <- rawGenomeFootprintMetrics[promoterIdx,]
      promoterGenomeMotifSignal <- (mean(promoterGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterGenomeFlankSignal <- (mean(promoterGenomeFootprintMetrics[,2], trim = 0.10) / 100)
      promoterGenomeBackgroundSignal <- (mean(promoterGenomeFootprintMetrics[,1], trim = 0.10) / 100)
      promoterGenome.log2Flank <- log2(promoterGenomeFlankSignal / promoterGenomeBackgroundSignal)
      promoterGenome.log2Depth <- log2(promoterGenomeMotifSignal / promoterGenomeFlankSignal)
      
      ## distal genome Sites
      distalGenomeFootprintMetrics <- rawGenomeFootprintMetrics[-promoterIdx,]
      distalGenomeMotifSignal <- (mean(distalGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalGenomeFlankSignal <- (mean(distalGenomeFootprintMetrics[,2], trim = 0.10) / 100)
      distalGenomeBackgroundSignal <- (mean(distalGenomeFootprintMetrics[,1], trim = 0.10) / 100)
      distalGenome.log2Flank <- log2(distalGenomeFlankSignal / distalGenomeBackgroundSignal)
      distalGenome.log2Depth <- log2(distalGenomeMotifSignal / distalGenomeFlankSignal)
      
      ## promoter Bound
      promoterBoundFootprintMetrics <- rawGenomeFootprintMetrics[promoterBoundIdx,]
      promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
      promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
      promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
      promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)
      
      ## promoter unbound
      promoterUnboundFootprintMetrics <- rawGenomeFootprintMetrics[promoterUnboundIdx,]
      promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
      promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
      promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
      promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)
      
      ## distal Bound
      distalBoundFootprintMetrics <- rawGenomeFootprintMetrics[distalBoundIdx,]
      distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
      distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
      distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
      distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)
      
      ## distal unbound
      distalUnboundFootprintMetrics <- rawGenomeFootprintMetrics[distalUnboundIdx,]
      distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
      distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
      distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
      distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)
      
      ## STORE THE DATA
      processedFootprintData$"promoterGenomeSites" <- promoterGenomeSites
      processedFootprintData$"promoterGenomeFootprintMetrics" <- promoterGenomeFootprintMetrics
      processedFootprintData$"promoterGenomeMotifSignal" <- promoterGenomeMotifSignal
      processedFootprintData$"promoterGenomeFlankSignal" <- promoterGenomeFlankSignal
      processedFootprintData$"promoterGenomeBackgroundSignal" <- promoterGenomeBackgroundSignal
      processedFootprintData$"promoterGenome.log2Flank" <- promoterGenome.log2Flank
      processedFootprintData$"promoterGenome.log2Depth" <- promoterGenome.log2Depth
      ##
      processedFootprintData$"distalGenomeSites" <- distalGenomeSites
      processedFootprintData$"distalGenomeFootprintMetrics" <- distalGenomeFootprintMetrics
      processedFootprintData$"distalGenomeMotifSignal" <- distalGenomeMotifSignal
      processedFootprintData$"distalGenomeFlankSignal" <- distalGenomeFlankSignal
      processedFootprintData$"distalGenomeBackgroundSignal" <- distalGenomeBackgroundSignal
      processedFootprintData$"distalGenome.log2Flank" <- distalGenome.log2Flank
      processedFootprintData$"distalGenome.log2Depth" <- distalGenome.log2Depth
      ##
      processedFootprintData$"promoterBoundSites" <- promoterBoundSites
      processedFootprintData$"promoterBoundFootprintMetrics" <- promoterBoundFootprintMetrics
      processedFootprintData$"promoterBoundMotifSignal" <- promoterBoundMotifSignal
      processedFootprintData$"promoterBoundFlankSignal" <- promoterBoundFlankSignal
      processedFootprintData$"promoterBoundBackgroundSignal" <- promoterBoundBackgroundSignal
      processedFootprintData$"promoterBound.log2Flank" <- promoterBound.log2Flank
      processedFootprintData$"promoterBound.log2Depth" <- promoterBound.log2Depth
      ##
      processedFootprintData$"promoterUnboundSites" <- promoterUnboundSites
      processedFootprintData$"promoterUnboundFootprintMetrics" <- promoterUnboundFootprintMetrics
      processedFootprintData$"promoterUnboundMotifSignal" <- promoterUnboundMotifSignal
      processedFootprintData$"promoterUnboundFlankSignal" <- promoterUnboundFlankSignal
      processedFootprintData$"promoterUnboundBackgroundSignal" <- promoterUnboundBackgroundSignal
      processedFootprintData$"promoterUnbound.log2Flank" <- promoterUnbound.log2Flank
      processedFootprintData$"promoterUnbound.log2Depth" <- promoterUnbound.log2Depth
      ##
      processedFootprintData$"distalBoundSites" <- distalBoundSites
      processedFootprintData$"distalBoundFootprintMetrics" <- distalBoundFootprintMetrics
      processedFootprintData$"distalBoundMotifSignal" <- distalBoundMotifSignal
      processedFootprintData$"distalBoundFlankSignal" <- distalBoundFlankSignal
      processedFootprintData$"distalBoundBackgroundSignal" <- distalBoundBackgroundSignal
      processedFootprintData$"distalBound.log2Flank" <- distalBound.log2Flank
      processedFootprintData$"distalBound.log2Depth" <- distalBound.log2Depth
      ##
      processedFootprintData$"distalUnboundSites" <- distalUnboundSites
      processedFootprintData$"distalUnboundFootprintMetrics" <- distalUnboundFootprintMetrics
      processedFootprintData$"distalUnboundMotifSignal" <- distalUnboundMotifSignal
      processedFootprintData$"distalUnboundFlankSignal" <- distalUnboundFlankSignal
      processedFootprintData$"distalUnboundBackgroundSignal" <- distalUnboundBackgroundSignal
      processedFootprintData$"distalUnbound.log2Flank" <- distalUnbound.log2Flank
      processedFootprintData$"distalUnbound.log2Depth" <- distalUnbound.log2Depth
      
      ## Save the data
      save(processedFootprintData, file = dataOutPath)
      
    } else {
      
      motifNames <- names(footprintData)
      ## The insertion matrices cannot be concatenated because they have different numbers of columns
      ## Initialize list objects here to store the insertion matrices separately (not concatenated)
      ## Pull the data for each individual motif and perform the bound/unbound split
      cat("Processing multiple motifs", "\n")
      ##
      for (z in 1:numMotifs){
        ## Pull the basic data
        cat("Pulling genome sites", "\n")
        com <- paste0("genomeSites", z, " <- footprintData[['", motifNames[z], "']][['genomeSites']]")
        eval(parse(text = com))
        cat("Pulling insertion matrix", "\n")
        com <- paste0("genomeInsertionMatrix", z, " <- footprintData[['", motifNames[z], "']][['insMatrix']]")
        eval(parse(text = com))
        cat("Pulling footprint metrics", "\n")
        com <- paste0("rawGenomeFootprintMetrics", z, " <- footprintData[['", motifNames[z], "']][['rawFootprintMetrics']]")
        eval(parse(text = com))
        ## Perform the bound/unbound split
        cat("Finding bound overlaps", "\n")
        com <- paste0("boundSiteOverlaps", z, " <- findOverlaps(footprintData[['", motifNames[z], "']][['parseData']][['bfSites']], footprintData[['", motifNames[z], "']][['genomeSites']])")
        eval(parse(text = com))
        cat("Setting bound index", "\n")
        com <- paste0("boundSiteIndex", z, " <- boundSiteOverlaps", z, "@to")
        eval(parse(text = com))
        ## Bound sites
        cat("Setting bound site data", "\n")
        com <- paste0("boundSites", z, " <- genomeSites", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("boundSitesInsertionMatrix", z, " <- genomeInsertionMatrix", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics", z, " <- rawGenomeFootprintMetrics", z, "[boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("numBoundSites", z, " <- length(boundSites", z, ")")
        eval(parse(text = com))
        ## Unbound sites
        cat("Setting unbound site data", "\n")
        com <- paste0("unboundSites", z, " <- genomeSites", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("unboundSitesInsertionMatrix", z, " <- genomeInsertionMatrix", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics", z, " <- rawGenomeFootprintMetrics", z, "[-boundSiteIndex", z, ",]")
        eval(parse(text = com))
        com <- paste0("numUnboundSites", z, " <- length(unboundSites", z, ")")
        eval(parse(text = com))
      } # end for (z in 1:numMotif)
      
      #### FIX THIS CODE LATER ####
      #### ADJUSTING VALUES TO PER BP ####
      for (m in 1:numMotifs){
        
        com <- paste0("motifWidth <- footprintData[['", motifNames[m], "']][['motifWidth']]")
        eval(parse(text = com))
        
        com <- paste0("boundSitesMetrics", m, "[,1] <- boundSitesMetrics", m, "[,1] / 100")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics", m, "[,2] <- boundSitesMetrics", m, "[,2] / 100")
        eval(parse(text = com))
        com <- paste0("boundSitesMetrics", m, "[,3] <- boundSitesMetrics", m, "[,3] / motifWidth")
        eval(parse(text = com))
        
        com <- paste0("rawGenomeFootprintMetrics", m, "[,1] <- rawGenomeFootprintMetrics", m, "[,1] / 100")
        eval(parse(text = com))
        com <- paste0("rawGenomeFootprintMetrics", m, "[,2] <- rawGenomeFootprintMetrics", m, "[,2] / 100")
        eval(parse(text = com))
        com <- paste0("rawGenomeFootprintMetrics", m, "[,3] <- rawGenomeFootprintMetrics", m, "[,3] / motifWidth")
        eval(parse(text = com))
        
        com <- paste0("unboundSitesMetrics", m, "[,1] <- unboundSitesMetrics", m, "[,1] / 100")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics", m, "[,2] <- unboundSitesMetrics", m, "[,2] / 100")
        eval(parse(text = com))
        com <- paste0("unboundSitesMetrics", m, "[,3] <- unboundSitesMetrics", m, "[,3] / motifWidth")
        eval(parse(text = com))
      }
      
      #### Perform the merging and deduplication ####
      for (b in 2:numMotifs){
        ## Find overlaps and generate selection indices
        com <- paste0("tempOverlapsGenome <- findOverlaps(genomeSites1, genomeSites", b,")")
        eval(parse(text = com))
        com <- paste0("tempOverlapsBound <- findOverlaps(boundSites1, boundSites", b,")")
        eval(parse(text = com))
        com <- paste0("tempOverlapsUnbound <- findOverlaps(unboundSites1, unboundSites", b,")")
        eval(parse(text = com))
        
        #### MERGE THE GENOME SITES (all sites) ####
        ## If no overlaps are present, can just directly merge the two Granges
        ## Otherwise, if some overlaps are present, merge the Granges, but omit overlapping sites from second group
        if (length(tempOverlapsGenome@from) == 0){
          com <- paste0("genomeSites1 <- c(genomeSites1, genomeSites", b, ")")
          eval(parse(text = com))
          com <- paste0("rawGenomeFootprintMetrics1 <- rbind(rawGenomeFootprintMetrics1, rawGenomeFootprintMetrics", b, ")")
          eval(parse(text = com))
        } else {
          mergeIdx <- tempOverlapsGenome@to
          com <- paste0("genomeSites1 <- c(genomeSites1, genomeSites", b, "[-mergeIdx])")
          eval(parse(text = com))
          com <- paste0("rawGenomeFootprintMetrics1 <- rbind(rawGenomeFootprintMetrics1, rawGenomeFootprintMetrics", b, "[-mergeIdx,])")
          eval(parse(text = com))
        } # end if (length(overlaps@from) == 0)
        
        #### MERGE THE BOUND SITES ####
        if (length(tempOverlapsBound@from) == 0){
          com <- paste0("boundSites1 <- c(boundSites1, boundSites", b, ")")
          eval(parse(text = com))
          com <- paste0("boundSitesMetrics1 <- rbind(boundSitesMetrics1, boundSitesMetrics", b, ")")
          eval(parse(text = com))
        } else {
          mergeIdx <- tempOverlapsBound@to
          com <- paste0("boundSites1 <- c(boundSites1, boundSites", b, "[-mergeIdx])")
          eval(parse(text = com))
          com <- paste0("boundSitesMetrics1 <- rbind(boundSitesMetrics1, boundSitesMetrics", b, "[-mergeIdx,])")
          eval(parse(text = com))
        } # end if (length(tempOverlapsBound@from) == 0)
        
        #### MERGE THE UNBOUND SITES ####
        if (length(tempOverlapsUnbound@from) == 0){
          com <- paste0("unboundSites1 <- c(unboundSites1, unboundSites", b, ")")
          eval(parse(text = com))
          com <- paste0("unboundSitesMetrics1 <- rbind(unboundSitesMetrics1, unboundSitesMetrics", b, ")")
          eval(parse(text = com))
        } else {
          mergeIdx <- tempOverlapsUnbound@to
          com <- paste0("unboundSites1 <- c(unboundSites1, unboundSites", b, "[-mergeIdx])")
          eval(parse(text = com))
          com <- paste0("unboundSitesMetrics1 <- rbind(unboundSitesMetrics1, unboundSitesMetrics", b, "[-mergeIdx,])")
          eval(parse(text = com))
        } # end if (length(tempOverlapsUnbound@from) == 0)
      } # end for (b in 2:numMotifs)
      
      ## Because the code is written to use the first motif to merge everything into,
      ## transfer the data at this stage to make it consistent with numMotifs == 1
      genomeSites <- genomeSites1
      rawGenomeFootprintMetrics <- rawGenomeFootprintMetrics1
      numGenomeSites <- length(genomeSites)
      ##
      boundSites <- boundSites1
      boundSitesMetrics <- boundSitesMetrics1
      numBoundSites <- length(boundSites)
      ##
      unboundSites <- unboundSites1
      unboundSitesMetrics <- unboundSitesMetrics1
      numUnboundSites <- length(unboundSites)
      
      #### Calculate footprint characteristics on merged data ####
      ## Calculate the 10% trimmed mean of all insertions in the motif sites
      genomeMotifSignal <- mean(rawGenomeFootprintMetrics[,3], trim = 0.10)
      boundMotifSignal <- mean(boundSitesMetrics[,3], trim = 0.10)
      unboundMotifSignal <- mean(unboundSitesMetrics[,3], trim = 0.10)
      
      ## Calculate the mean of all insertions in the flank region
      genomeFlankSignal <- mean(rawGenomeFootprintMetrics[,2], trim = 0.10)
      boundFlankSignal <- mean(boundSitesMetrics[,2], trim = 0.10)
      unboundFlankSignal <- mean(unboundSitesMetrics[,2], trim = 0.10)
      
      ## Calculate the mean of background insertions
      genomeBackgroundSignal <- mean(rawGenomeFootprintMetrics[,1], trim = 0.10)
      boundBackgroundSignal <- mean(boundSitesMetrics[,1], trim = 0.10)
      unboundBackgroundSignal <- mean(unboundSitesMetrics[,1], trim = 0.10)
      
      ## Calculate flanking accessibility (log2 fold change between flank and background)
      genome.log2Flank <- log2(genomeFlankSignal / genomeBackgroundSignal)
      bound.log2Flank <- log2(boundFlankSignal / boundBackgroundSignal)
      unbound.log2Flank <- log2(unboundFlankSignal / unboundBackgroundSignal)
      
      ## Calculate footprint depth (log2 fold change between flank and background)
      genome.log2Depth <- log2(genomeMotifSignal / genomeFlankSignal)
      bound.log2Depth <- log2(boundMotifSignal / boundFlankSignal)
      unbound.log2Depth <- log2(unboundMotifSignal / unboundFlankSignal)
      
      #### TRANSFER DATA TO STORAGE OBJECT ####
      ## Initialize a new list object to store the processed data
      processedFootprintData <- list()
      ##
      processedFootprintData$"geneName" <- geneName
      processedFootprintData$"numMotifs" <- numMotifs
      processedFootprintData$"numGenomeSites" <- numGenomeSites
      processedFootprintData$"numBoundSites" <- numBoundSites
      processedFootprintData$"numUnboundSites" <- numUnboundSites
      ##
      processedFootprintData$"genomeSites" <- genomeSites
      processedFootprintData$"rawGenomeFootprintMetrics" <- rawGenomeFootprintMetrics
      processedFootprintData$"genomeMotifSignal" <- genomeMotifSignal
      processedFootprintData$"genomeFlankSignal" <- genomeFlankSignal
      processedFootprintData$"genomeBackgroundSignal" <- genomeBackgroundSignal
      processedFootprintData$"genome.log2Flank" <- genome.log2Flank
      processedFootprintData$"genome.log2Depth" <- genome.log2Depth
      ##
      processedFootprintData$"boundSites" <- boundSites
      processedFootprintData$"boundSitesMetrics" <- boundSitesMetrics
      processedFootprintData$"boundMotifSignal" <- boundMotifSignal
      processedFootprintData$"boundFlankSignal" <- boundFlankSignal
      processedFootprintData$"boundBackgroundSignal" <- boundBackgroundSignal
      processedFootprintData$"bound.log2Flank" <- bound.log2Flank
      processedFootprintData$"bound.log2Depth" <- bound.log2Depth
      ##
      processedFootprintData$"unboundSites" <- unboundSites
      processedFootprintData$"unboundSitesMetrics" <- unboundSitesMetrics
      processedFootprintData$"unboundMotifSignal" <- unboundMotifSignal
      processedFootprintData$"unboundFlankSignal" <- unboundFlankSignal
      processedFootprintData$"unboundBackgroundSignal" <- unboundBackgroundSignal
      processedFootprintData$"unbound.log2Flank" <- unbound.log2Flank
      processedFootprintData$"unbound.log2Depth" <- unbound.log2Depth
      
      #### CODE TESTING - PROMOTER/DISTAL groups ####
      ## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
      promoters <- promoters(txdb, upstream = 1000, downstream = 100)
      ## Trim the GRanges object to keep standard entries only
      scope <- paste0("chr", c(1:22, "X", "Y"))
      promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
      promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
      
      ## Subset based on the overlaps
      promoterOverlaps <- findOverlaps(promoters, genomeSites, ignore.strand = TRUE)
      promoterIdx <- unique(promoterOverlaps@to)
      ##
      promoterGenomeSites <- genomeSites[promoterIdx]
      distalGenomeSites <- genomeSites[-promoterIdx]
      ##
      promoterBoundOverlaps <- findOverlaps(promoterGenomeSites, boundSites)
      promoterUnboundOverlaps <- findOverlaps(promoterGenomeSites, unboundSites)
      promoterBoundIdx <- unique(promoterBoundOverlaps@from)
      promoterUnboundIdx <- unique(promoterUnboundOverlaps@from)
      ##
      promoterBoundSites <- genomeSites[promoterBoundIdx]
      promoterUnboundSites <- genomeSites[promoterUnboundIdx]
      ##
      distalBoundOverlaps <- findOverlaps(distalGenomeSites, boundSites)
      distalUnboundOverlaps <- findOverlaps(distalGenomeSites, unboundSites)
      distalBoundIdx <- unique(distalBoundOverlaps@from)
      distalUnboundIdx <- unique(distalUnboundOverlaps@from)
      ##
      distalBoundSites <- genomeSites[distalBoundIdx]
      distalUnboundSites <- genomeSites[distalUnboundIdx]
      
      ##### promoter genome sites
      promoterGenomeFootprintMetrics <- rawGenomeFootprintMetrics[promoterIdx,]
      promoterGenomeMotifSignal <- (mean(promoterGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterGenomeFlankSignal <- (mean(promoterGenomeFootprintMetrics[,2], trim = 0.10) / 100)
      promoterGenomeBackgroundSignal <- (mean(promoterGenomeFootprintMetrics[,1], trim = 0.10) / 100)
      promoterGenome.log2Flank <- log2(promoterGenomeFlankSignal / promoterGenomeBackgroundSignal)
      promoterGenome.log2Depth <- log2(promoterGenomeMotifSignal / promoterGenomeFlankSignal)
      
      ## distal genome Sites
      distalGenomeFootprintMetrics <- rawGenomeFootprintMetrics[-promoterIdx,]
      distalGenomeMotifSignal <- (mean(distalGenomeFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalGenomeFlankSignal <- (mean(distalGenomeFootprintMetrics[,2], trim = 0.10) / 100)
      distalGenomeBackgroundSignal <- (mean(distalGenomeFootprintMetrics[,1], trim = 0.10) / 100)
      distalGenome.log2Flank <- log2(distalGenomeFlankSignal / distalGenomeBackgroundSignal)
      distalGenome.log2Depth <- log2(distalGenomeMotifSignal / distalGenomeFlankSignal)
      
      ## promoter Bound
      promoterBoundFootprintMetrics <- rawGenomeFootprintMetrics[promoterBoundIdx,]
      promoterBoundMotifSignal <- (mean(promoterBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterBoundFlankSignal <- (mean(promoterBoundFootprintMetrics[,2], trim = 0.10) / 100)
      promoterBoundBackgroundSignal <- (mean(promoterBoundFootprintMetrics[,1], trim = 0.10) / 100)
      promoterBound.log2Flank <- log2(promoterBoundFlankSignal / promoterBoundBackgroundSignal)
      promoterBound.log2Depth <- log2(promoterBoundMotifSignal / promoterBoundFlankSignal)
      
      ## promoter unbound
      promoterUnboundFootprintMetrics <- rawGenomeFootprintMetrics[promoterUnboundIdx,]
      promoterUnboundMotifSignal <- (mean(promoterUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      promoterUnboundFlankSignal <- (mean(promoterUnboundFootprintMetrics[,2], trim = 0.10) / 100)
      promoterUnboundBackgroundSignal <- (mean(promoterUnboundFootprintMetrics[,1], trim = 0.10) / 100)
      promoterUnbound.log2Flank <- log2(promoterUnboundFlankSignal / promoterUnboundBackgroundSignal)
      promoterUnbound.log2Depth <- log2(promoterUnboundMotifSignal / promoterUnboundFlankSignal)
      
      ## distal Bound
      distalBoundFootprintMetrics <- rawGenomeFootprintMetrics[distalBoundIdx,]
      distalBoundMotifSignal <- (mean(distalBoundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalBoundFlankSignal <- (mean(distalBoundFootprintMetrics[,2], trim = 0.10) / 100)
      distalBoundBackgroundSignal <- (mean(distalBoundFootprintMetrics[,1], trim = 0.10) / 100)
      distalBound.log2Flank <- log2(distalBoundFlankSignal / distalBoundBackgroundSignal)
      distalBound.log2Depth <- log2(distalBoundMotifSignal / distalBoundFlankSignal)
      
      ## distal unbound
      distalUnboundFootprintMetrics <- rawGenomeFootprintMetrics[distalUnboundIdx,]
      distalUnboundMotifSignal <- (mean(distalUnboundFootprintMetrics[,3], trim = 0.10) / motifWidth)
      distalUnboundFlankSignal <- (mean(distalUnboundFootprintMetrics[,2], trim = 0.10) / 100)
      distalUnboundBackgroundSignal <- (mean(distalUnboundFootprintMetrics[,1], trim = 0.10) / 100)
      distalUnbound.log2Flank <- log2(distalUnboundFlankSignal / distalUnboundBackgroundSignal)
      distalUnbound.log2Depth <- log2(distalUnboundMotifSignal / distalUnboundFlankSignal)
      
      ## STORE THE DATA
      processedFootprintData$"promoterGenomeSites" <- promoterGenomeSites
      processedFootprintData$"promoterGenomeFootprintMetrics" <- promoterGenomeFootprintMetrics
      processedFootprintData$"promoterGenomeMotifSignal" <- promoterGenomeMotifSignal
      processedFootprintData$"promoterGenomeFlankSignal" <- promoterGenomeFlankSignal
      processedFootprintData$"promoterGenomeBackgroundSignal" <- promoterGenomeBackgroundSignal
      processedFootprintData$"promoterGenome.log2Flank" <- promoterGenome.log2Flank
      processedFootprintData$"promoterGenome.log2Depth" <- promoterGenome.log2Depth
      ##
      processedFootprintData$"distalGenomeSites" <- distalGenomeSites
      processedFootprintData$"distalGenomeFootprintMetrics" <- distalGenomeFootprintMetrics
      processedFootprintData$"distalGenomeMotifSignal" <- distalGenomeMotifSignal
      processedFootprintData$"distalGenomeFlankSignal" <- distalGenomeFlankSignal
      processedFootprintData$"distalGenomeBackgroundSignal" <- distalGenomeBackgroundSignal
      processedFootprintData$"distalGenome.log2Flank" <- distalGenome.log2Flank
      processedFootprintData$"distalGenome.log2Depth" <- distalGenome.log2Depth
      ##
      processedFootprintData$"promoterBoundSites" <- promoterBoundSites
      processedFootprintData$"promoterBoundFootprintMetrics" <- promoterBoundFootprintMetrics
      processedFootprintData$"promoterBoundMotifSignal" <- promoterBoundMotifSignal
      processedFootprintData$"promoterBoundFlankSignal" <- promoterBoundFlankSignal
      processedFootprintData$"promoterBoundBackgroundSignal" <- promoterBoundBackgroundSignal
      processedFootprintData$"promoterBound.log2Flank" <- promoterBound.log2Flank
      processedFootprintData$"promoterBound.log2Depth" <- promoterBound.log2Depth
      ##
      processedFootprintData$"promoterUnboundSites" <- promoterUnboundSites
      processedFootprintData$"promoterUnboundFootprintMetrics" <- promoterUnboundFootprintMetrics
      processedFootprintData$"promoterUnboundMotifSignal" <- promoterUnboundMotifSignal
      processedFootprintData$"promoterUnboundFlankSignal" <- promoterUnboundFlankSignal
      processedFootprintData$"promoterUnboundBackgroundSignal" <- promoterUnboundBackgroundSignal
      processedFootprintData$"promoterUnbound.log2Flank" <- promoterUnbound.log2Flank
      processedFootprintData$"promoterUnbound.log2Depth" <- promoterUnbound.log2Depth
      ##
      processedFootprintData$"distalBoundSites" <- distalBoundSites
      processedFootprintData$"distalBoundFootprintMetrics" <- distalBoundFootprintMetrics
      processedFootprintData$"distalBoundMotifSignal" <- distalBoundMotifSignal
      processedFootprintData$"distalBoundFlankSignal" <- distalBoundFlankSignal
      processedFootprintData$"distalBoundBackgroundSignal" <- distalBoundBackgroundSignal
      processedFootprintData$"distalBound.log2Flank" <- distalBound.log2Flank
      processedFootprintData$"distalBound.log2Depth" <- distalBound.log2Depth
      ##
      processedFootprintData$"distalUnboundSites" <- distalUnboundSites
      processedFootprintData$"distalUnboundFootprintMetrics" <- distalUnboundFootprintMetrics
      processedFootprintData$"distalUnboundMotifSignal" <- distalUnboundMotifSignal
      processedFootprintData$"distalUnboundFlankSignal" <- distalUnboundFlankSignal
      processedFootprintData$"distalUnboundBackgroundSignal" <- distalUnboundBackgroundSignal
      processedFootprintData$"distalUnbound.log2Flank" <- distalUnbound.log2Flank
      processedFootprintData$"distalUnbound.log2Depth" <- distalUnbound.log2Depth
      
      ## Save the data
      save(processedFootprintData, file = dataOutPath)
      
    } # end if (numMotifs == 0)
  } # end if (file.exists(dataOutPath) == TRUE)
  
}, # end try
error=function(cond){
  message(cond)
  return(NA)
},
finally={})

##
file.create(outPath)
cat("Finished!", "\n")


#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

#### Library loading
cat("Loading libraries...", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))

##
cat("Setting snakemake vars...", "\n")
motifDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
rdataPath <- gsub("operations/scans", "data", outPath)
rdataPath <- gsub("PWMscan.done", "bindingSites.Rdata", rdataPath)
currentGene <- snakemake@wildcards[["gene"]]

##
cat("Loading motifData object...", "\n")
load(motifDataPath)

##
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))

## Set parameters
numMotifs <- length(motifs)
genome <- Hsapiens
score <- "99%"
bindingSites <- list()
cat("Found", numMotifs, " unique motifs for gene", currentGene, "scanning for sites with", score, "PWM match", "\n")

## Scan the genome for matches to each unique motif
for (a in 1:numMotifs){
  cat("Scanning motif", a, "for gene", currentGene, "\n")
  tempSites <- list()
  PWM <- motifs[[a]]
  sites <- matchPWM(PWM, genome, min.score=score)
  tempSites$PWM <- PWM
  tempSites$sites <- sites
  bindingSites[[a]] <- tempSites
} # end for (a in 1:numMotifs)

## Save the data
cat("Saving data...", "\n")
save(bindingSites, file = rdataPath)
file.create(outPath)

rm(list=ls())
gc()
cat("Finished scanning!", "\n")

##
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(MotifDb))

##
cat("Setting snakemake vars...", "\n")
sample1Data <- snakemake@input[[1]]
sample2Data <- snakemake@input[[2]]
outputPath <- snakemake@output[[1]]
sampleName1 <- snakemake@wildcards[["xsample1"]]
sampleName2 <- snakemake@wildcards[["xsample2"]]
geneName <- snakemake@wildcards[["gene"]]
motifNum <- snakemake@wildcards[["motif"]]

##
cat("Analyzing samples", sampleName1, sampleName2, "for gene", geneName, "motif number", motifNum, "\n")

## Set PWM for this motif
PWMin <- pwm2pfm(PWM)

## Load in the data from both samples
cat("Loading data for first sample from path:", sample1Data, "\n")
load(sample1Data)
cat("Loading data for second sample from path:", sample2Data, "\n")
load(sample2Data)

## Sample 1 setup
sample1Signals <-
  sample1Sites <- 
  
  ## Sample 2 setup
  sample2Signals <-
  sample2Sites <- 
  
  ## Find GRanges intersection for bf peak passing sites for both samples
  xSitesOverlap <- findOverlaps(sample1Sites, sample2Sites)
#xSites <- GenomicRanges::intersect(sample1Sites, sample2Sites)

## Subset the signals and sites objects for both samples based on the intersection
idx1 <- xSitesOverlap@from
idx2 <- xSitesOverlap@to


## Get common sites
sample1xSites <- sample1Sites[idx1]
sample2xSites <- sample2Sites[idx2]


## Get common signals
sample1xSignals <- sample1Signals[idx1,]
sample2xSignals <- sample2Signals[idx2,]




## Make sure the data is in the same order for both samples



## Generate and save footprint plots
svgPlot1 <- paste0("xsample/", sampleName1, ".", sampleName2, "/footprints/", sampleName1, ".", geneName, ".", motifNum, ".bfpeak.xsites.svg")
cat("Saving footprint plot for sample 1 at path:", svgPlot1, "\n")

##
svg(file = svgPlot1) # set the filepath for saving the svg figure
ATACseqQC:::plotFootprints(prof,
                           Mlen=motifWidth, motif=PWMin)
dev.off()



## Generate and save heatmap footprint plots

cat("Generating heatmap...", "\n")
cat("1", "\n")
heatsigs <- combinedbfPassPeakSignal
cat("2", "\n")
heatnumsites <- numbfPassPeakSites
cat("3", "\n")
heatnumbp <- length(combinedbfPassPeakSignal[1,])
cat("4", "\n")
heatsites <- bfPassPeakSites


#sigs <- parsedSitesInfo[["bfPassPeakSignals"]]
#sig_plus <- sigs["+"]
#sig_minus <- sigs["-"]
#numsites <- length(sig_plus[["+"]][,1])
#numbp <- length(sig_plus[["+"]][1,])
#sites <- parsedSitesInfo[["bfPassPeakSites"]]

##
#temp <- matrix(data = NA, nrow = numsites, ncol = numbp)
#for (i in 1:numsites){temp[i,] <- sigs[['+']][i,] + sigs[['-']][i,]}

## scale each row individually
cat("5", "\n")
for (a in 1:heatnumsites){
  maxsig <- max(heatsigs[a,])
  for (b in 1:heatnumbp){heatsigs[a,b] <- (heatsigs[a,b] / maxsig)}}
maxsig <- 1

## invert signals
cat("6", "\n")
for (a in 1:heatnumsites){for (b in 1:heatnumbp){heatsigs[a,b] <- (1-heatsigs[a,b])}}

## Annotate the combined sublist name which will become the title of the heatmap plot
cat("7", "\n")
plottitle <- paste0(genename, "_motif", x, "_numsites", heatnumsites)
combined <- list()
com <- paste0("combined$", plottitle, " <- heatsigs")
eval(parse(text = com))

##
cat("8", "\n")
heatsvg_path <- paste0(dirpath, "footprints/heatmaps/", samplename, ".", genename, ".", "motif", x, ".bfpeak.sites.heatmap.svg")
svg(file = heatsvg_path) # set the filepath for saving the svg figure
cat("Saving svg footprint image at path:", heatsvg_path, "\n")

## Margin controls
# margin(a,b,c,d)
# a = size of graph from top to bottom, higher value = smaller. default = 0.1
# b = size of graph from left to right, higher value = smaller. default = 0.005
# c = flips x axis?
# d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
# good settings for ATACseq = c(0.1,0.005,0.05,0.2)
# bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values

##
cat("9", "\n")
ChIPpeakAnno::featureAlignedHeatmap(combined,
                                    feature.gr=reCenterPeaks(heatsites,width=heatnumbp),
                                    upper.extreme = maxsig, # set this to control the heatmap scale
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=heatnumbp,
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias=0.9)(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()
############################################################











