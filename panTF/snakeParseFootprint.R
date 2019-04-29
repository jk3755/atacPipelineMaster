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
  #cat("Generating a null footprint model with the following parameters:", "\n")
  #cat("Iterations:", iterations, "\n")
  #cat("Input signal:", inputSignal, "\n")
  #cat("Analysis window (bp):", analysisWidth, "\n")
  #cat("Motif width (bp):", motifWidth, "\n")
  
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
    null <- c(as.vector(rmultinom(1, size=inputSignal, prob=rep(1, length(null)))))
    
    ## Calculate the mean signal in motif region
    motifStart <- ((analysisWidth - motifWidth)/2)
    motifEnd <- (motifStart + motifWidth)
    motifAvg <- (sum(null[motifStart:motifEnd])) / motifWidth
    
    ## Store the average values
    averages[a] <- motifAvg
    
  } # end for (a in 1:n)
  return(averages)
} # end generateNullFP function

## Set snakemake variables
cat("Setting snakemake variables", "\n")
footprintDataPath <- snakemake@input[[1]]
peaksPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Set the output filepath for the Rdata object and perform a filecheck
dataOutPath <- gsub("operations", "data", outPath)
dataOutPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", dataOutPath, perl = TRUE)
cat("Output path for parsed footprint data file set as:", dataOutPath, "\n")

## Perform a filecheck, and skip if file already exists
if (file.exists(dataOutPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  ## Load libraries
  cat("Loading libraries", "\n")
  suppressMessages(library(BiocGenerics))
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(GenomicAlignments))
  suppressMessages(library(stats4))
  suppressMessages(library(Rsamtools))
  suppressMessages(library(genomation))
  suppressMessages(library(seqLogo))
  suppressMessages(library(ChIPpeakAnno))
  suppressMessages(library(rlist))
  
  ## Load the footprintData file
  footprintDataPath <- gsub("operations", "data", footprintDataPath)
  footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
  cat("Loading footprintData file at path:", footprintDataPath, "\n")
  load(footprintDataPath)
  
  ## To avoid errors, clear the list of any empty sub-lists first
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  
  ## If the data object is empty, skip the parse operation and output a dummy file
  if (length(footprintData) == 0){
    cat("No data found in footprint object. Skipping", "\n")
  } else {
    
    ## Performing parsing operations
    numMotif <- length(footprintData)
    motifNames <- names(footprintData)
    cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
    
    ## Load the peaks file
    cat("Loading accessibility peaks and trimming ranges", "\n")
    scope <- paste0("chr", c(1:22, "X", "Y"))
    grPeaks <- readBed(peaksPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
    grPeaks <- keepStandardChromosomes(grPeaks, pruning.mode="coarse")
    grPeaks <- keepSeqlevels(grPeaks, scope, pruning.mode="coarse")
    grPeaks <- trim(grPeaks, use.names = TRUE)
    
    ## Loop over all motifs with data
    for (a in 1:numMotif){
      
      ##
      cat("Processing motif", a, "\n")
      #tryCatch({
        
        ############################################################
        ## Load the temporary data
        com <- paste0("tempData <- footprintData$", motifNames[a])
        eval(parse(text = com))
        ##
        librarySize <- tempData[["librarySize"]]
        coverageSize <- tempData[["coverageSize"]]
        libraryFactor <- tempData[["libraryFactor"]]
        PWM <- tempData[["PWM"]]
        genomeSites <- tempData[["genomeSites"]]
        numGenomeSites <- tempData[["numGenomeSites"]]
        motifWidth <- tempData[["motifWidth"]]
        extendedSites <- tempData[["extendedSites"]]
        shiftedInsertions <- tempData[["shiftedInsertions"]]
        insertionRLE <- tempData[["insertionRLE"]]
        insertionViews <- tempData[["insertionViews"]]
        insertionMatrix <- tempData[["insertionMatrix"]]
        ############################################################
        
        ## Process all the binding sites, use null model
        ## Determine which sites are bound and which are not
        cat("Processing binding sites to find bound TFs", "\n")
        
        
        insVector <- insProfile[2,]
        siteTotalSignal <- c()
        
        
        
        ## Calculate total signal for each site
        for (b in 1:numSites){
          siteTotalSignal[b] <- sum(tempData[["insMatrix"]][b,])} # end for (b in 1:numSites)
        
        ## Find the unique values for total signal and generate null models
        uniqueTotalSignals <- unique(siteTotalSignal)
        ## Remove NA values from uniqueTotalSignals
        ## (how do they get there???)
        uniqueTotalSignals <- uniqueTotalSignals[!is.na(uniqueTotalSignals)]
        
        ## Initiate a matrix to store the mean null signal in the null model and the input signal to null model
        nullModels <- matrix(data = NA, ncol = 2, nrow = length(uniqueTotalSignals))
        colnames(nullModels) <- c("Input signal", "Avg motif signal in null model")
        
        ## Calculate the null models
        for (c in 1:length(uniqueTotalSignals)){
          nullVec <- generateNullFP(1000, uniqueTotalSignals[c], siteWidth, motifWidth)
          nullModels[c,1] <- uniqueTotalSignals[c]
          nullModels[c,2] <- mean(nullVec)} # end for (c in 1:length(uniqueTotalSignals))
        
        ## Perform a one-tailed t-test to generate a p-value for each observed motif site
        cat("Performing one-tailed t-tests on peak subset", "\n")
        ttestPeak <- list() # list to store the results of the t-tests
        pvaluePeak <- c() # vector to store the p-values
        tvaluePeak <- c() # vector to store the t-value
        
        ## Perform t-test on all sites
        for (d in 1:numSites){
          ## Retrieve the total signal for the current site
          currentSignal <- c(siteTotalSignal[d])
          ## Retrieve the appropriate null model
          currentNullModel <- nullModels[which(nullModels[,1]==currentSignal),2]
          ## Perform the t-test
          ttestPeak[[d]] <- t.test(insMatrix[d,250:(250+motifWidth)], mu=currentNullModel, alternative="less", conf.level = 0.95)
          pvaluePeak[d] <- ttestPeak[[d]][["p.value"]]
          tvaluePeak[d] <- ttestPeak[[d]][["statistic"]][["t"]]
        } # for (d in 1:numSites)
        
        ## Get the indices of the sites that are lower than p = 0.05
        cat("Selecting p-value passing sites", "\n")
        idxPvaluePass <- which(pvaluePeak < 0.05)
        peakPvaluePass <- pvaluePeak[idxPvaluePass]
        
        ## Perform bonferroni correction
        cat("Performing bonferroni correction", "\n")
        idxbfPeakPass <- which(pvaluePeak < (0.05/numSites))
        bfPvaluePeakPass <- pvaluePeak[idxbfPeakPass]
        
        ## Subset the insertion matrix based on the bonferroni passing sites only
        cat("Subsetting sites based on bf corrected p-values", "\n")
        bfInsMatrix <- insMatrix[idxbfPeakPass,]
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
        bfSites <- peakSites[idxbfPeakPass]
        bfNumSites <- length(idxbfPeakPass)
        
        ## Calculate the insertion probability at each basepair
        cat("Calculating insertion probabilies", "\n")
        rawTotalSignal <- sum(insertionMatrix)
        rawProfile <- matrix(data = NA, ncol = length(insertionMatrix[1,]), nrow = 2)
        rownames(rawProfile) <- c("Column sums", "Insertion frequency")
        ##
        for (c in 1:length(insertionMatrix[1,])){
          rawProfile[1,c] <- sum(insertionMatrix[,c])
          rawProfile[2,c] <- (rawProfile[1,c] / rawTotalSignal) * 100
        } # end for (c in 1:length(insMatrix[1,]))
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
        ## To avoid errors, clear the list of any empty sub-lists first
        ## Should this result in an object with no data, that can be output
        ## as a dummy file to keep the pipeline running smoothly
        footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
        ### ADD THIS CODE HERE INSTEAD OF IN RAW ANALYSIS CODE ########################################################
        
        
        
        
        
        
        
        ## Analyze the binding sites in peaks
        cat("Processing binding sites in accessibility peaks", "\n")
        peakOverlaps <- findOverlaps(genomeSites, grPeaks)
        peakIndex <- unique(peakOverlapsIdx@to)
        peakSites <- genomeSites[peakIndex]
        peakSites <- keepStandardChromosomes(peakSites, pruning.mode="coarse")
        peakSites <- keepSeqlevels(peakSites, scope, pruning.mode="coarse")
        peakSites <- trim(peakSites, use.names = TRUE)
        numPeakSites <- length(peakSites)
        cat("Found", numPeakSites, "motif binding sites in peak accessibility regions", "\n")
        ##
        peakInsertionMatrix <- insertionMatrix[peakIndex,]
        
        
        
        
        
        ## Data transfer to storage object and save
        cat("Transferring data to storage object", "\n")
        parseData <- list()
        ##
        parseData$numSites <- numSites
        parseData$insVector <- insVector
        parseData$siteTotalSignal <- siteTotalSignal
        parseData$uniqueTotalSignals <- uniqueTotalSignals
        parseData$nullModels <- nullModels
        parseData$ttestPeak <- ttestPeak
        parseData$pvaluePeak <- pvaluePeak
        parseData$tvaluePeak <- tvaluePeak
        parseData$peakPvaluePass <- peakPvaluePass
        parseData$bfPvaluePeakPass <- bfPvaluePeakPass
        parseData$bfInsMatrix <- bfInsMatrix
        parseData$bfTotalSignal <- bfTotalSignal
        parseData$bfProfile <- bfProfile
        parseData$bfVector <- bfVector
        parseData$bfSites <- bfSites
        parseData$bfNumSites <- bfNumSites
        ##
        com <- paste0("footprintData$motif", a, "$parseData <- parseData")
        eval(parse(text = com))
        
      # }, # end try
      # error=function(cond){
      #   message(cond)
      #   return(NA)
      # },
      # finally={})

    } # end for (a in 1:numMotif)
  } # end if (length(footprintData) == 0)
  
  ## Save the data
  save(footprintData, file = dataOutPath)
  
} # end if (file.exists(dataOutPath) == TRUE)






## Create the output file for snakemake
file.create(outPath)
cat("Finished!", "\n")



