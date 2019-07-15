## The processing script should do the following things:
# Parse the footprint sites into peaks / non-peaks
# Use Null model to determine binding or not binding at every individual site genome-wide
# 

## Disable scientific notation in variables
options(scipen = 999)
## suppress warnings globally here, as they will disrupt the tryCatch block
## will need to improve this code at some point
options(warn = -1)

#### Set snakemake variables
footprintDataPath <- snakemake@input[[1]]
sampleTotalReadsPath <- snakemake@input[[2]]
peaksPath <- snakemake@input[[3]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

#### Set hg38 number of bases (haploid)
hg38TotalBP <- 3272116950

#### Set the output filepath for the Rdata object and perform a filecheck
cat("Processing footprint data for", geneName, "from sample", sampleName, "\n")
dataOutPath <- gsub("operations/parse", "data/parsed", outPath)
dataOutPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", dataOutPath, perl = TRUE)
cat("Output path for parsed data:", dataOutPath, "\n")

#### Perform a filecheck 
if (file.exists(dataOutPath) == TRUE){
  cat("File already exists, skipping", "\n")
} else {
  
  #### Load libraries
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
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  #### Load hg38 annotations
  cat("Loading hg38 annotations from TxDb", "\n")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb <- keepStandardChromosomes(txdb, pruning.mode="coarse")
  
  #### Build functions
  cat("Building functions...", "\n")
  generateNullFP <- function(iterations, inputSignal, analysisWidth, motifWidth){
    # This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
    # Consider the total signal (number of insertions) at each specific ~200 bp locus
    # Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
    # Use published or experimentally derived models of Tn5 sequence specific insertion bias
    # For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
    # Generate the null model graph by weighted random residstribution of the total observed signal at that site
    # Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
    # These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
    # iterations = number of iterations
    # inputSignals = unique values for total signal
    # analysisWidth = total bp in region of interest (flank + background + motif)
    # motifWidth = motif width
    
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
  
  #### Load the raw footprintData file
  footprintDataPath <- gsub("operations", "data", footprintDataPath)
  footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
  cat("Loading raw footprintData file from path:", footprintDataPath, "\n")
  load(footprintDataPath)
  
  #### To avoid errors, clear the list of any empty sub-lists first
  ## Remove motifs with 0 binding sites
  cat("Removing motifs that matched 0 genomic loci", "\n")
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  numMotif <- length(footprintData)
  cat("Found data for", numMotif, "motifs", "\n")
  
  #### If no raw data is found, skip this footprint
  if (numMotif == 0){
    cat(numMotif, "No data to analyze. Skipping", "\n")
  } else {
    
    #### Get the total number of reads in the sample
    cat("Loading total sample reads from:", sampleTotalReadsPath, "\n")
    load(sampleTotalReadsPath)
    cat("Found", sampleTotalReads, "total reads in current sample", "\n")
    
    #### Load the peaks data for current sample
    cat("Loading sample accesibility peak data from:", peaksPath, "\n")
    samplePeaks <- readBed(peaksPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
    samplePeaks <- keepStandardChromosomes(samplePeaks, pruning.mode="coarse")
    samplePeaks <- trim(samplePeaks)
    
    #### Loop over all motifs, perform the processing operations
    cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
    for (a in 1:numMotif){
      
      #### Pull the data from the raw footprint analysis
      cat("Processing footprint data for motif", a, "\n")
      com <- paste0("tempData <- footprintData$motif", a)
      eval(parse(text = com))
      ##
      PWM <- tempData[["PWM"]]
      motifWidth <- tempData[["motifWidth"]]
      allSites <- tempData[["allSites"]]
      extSites <- tempData[["extSites"]]
      insMatrix <- tempData[["insMatrix"]]
      libSize <- tempData[["libsize"]]
      coverageSize <- tempData[["coverageSize"]]
      libFactor <- tempData[["libFactor"]]
      ## Note that because trimming to standard xsomes is done, will need to set
      ## the total number of sites to the row number of the insertion matrix here
      numSites <- length(insMatrix[,1])
      
      #### Calculate basic statistics for each site with raw data
      rawSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 10)
      colnames(rawSiteBasicStats) <- c("Site index", "Total signal", "Total signal per bp", "Motif signal per bp",
                                    "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                    "Flank / Background", "Motif / Flank", "Motif / Wide Flank") 
      
      #### Populate the basic stats matrix
      for (b in 1:numSites){
        rawSiteBasicStats[b,1] <- b # Site index
        rawSiteBasicStats[b,2] <- sum(insMatrix[b,]) # Total signal
        rawSiteBasicStats[b,3] <- rawSiteBasicStats[b,2] / (500 + motifWidth) # Total signal per bp
        rawSiteBasicStats[b,4] <- sum(insMatrix[b,(250:(250 + motifWidth))] / motifWidth) # Motif signal per bp
        rawSiteBasicStats[b,5] <- (sum(insMatrix[b,200:250]) + sum(insMatrix[b,(250 + motifWidth):(300 + motifWidth)])) / 100 # Flank signal per bp
        rawSiteBasicStats[b,6] <- (sum(insMatrix[b,1:50]) + sum(insMatrix[b,(500 + motifWidth-50):(500 + motifWidth)])) / 100 # Background signal per bp
        rawSiteBasicStats[b,7] <- (sum(insMatrix[b,1:250]) + sum(insMatrix[b,(250 + motifWidth):(500 + motifWidth)])) / 500 # Wide flank signal per bp
        rawSiteBasicStats[b,8] <- rawSiteBasicStats[b,5] / rawSiteBasicStats[b,6] # Flank / background
        rawSiteBasicStats[b,9] <- rawSiteBasicStats[b,4] / rawSiteBasicStats[b,5] # Motif / flank
        rawSiteBasicStats[b,10] <- rawSiteBasicStats[b,4] / rawSiteBasicStats[b,7] # Motif / wide flank
      } # end for (b in 1:numSites)
      
      #### Generate the insertion site probability vector for raw data
      rawInsProb <- c()
      for (c in 1:(500 + motifWidth)){
        rawInsProb[c] <- sum(insMatrix[,c])
      } # end for (c in 1:(500 + motifWidth))
      rawTotalSignal<- sum(rawInsProb)
      rawInsProb <- rawInsProb / rawTotalSignal
      
      #### Make the CPM normalized insertion matrix ####
      ## This is an important normalization that will allow you to compare one sample to another directly,
      ## Even if the two samples have a different number of total reads
      
      ##
      factorCPM <- sampleTotalReads / 1000000
      CPMNinsMatrix <- insMatrix / factorCPM
      
      ##
      #### Calculate basic statistics for each site with raw data
      CPMNSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 10)
      colnames(CPMNSiteBasicStats) <- c("Site index", "Total signal", "Total signal per bp", "Motif signal per bp",
                                       "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                       "Flank / Background", "Motif / Flank", "Motif / Wide Flank") 
      
      #### Populate the basic stats matrix
      for (b in 1:numSites){
        CPMNSiteBasicStats[b,1] <- b # Site index
        CPMNSiteBasicStats[b,2] <- sum(CPMNinsMatrix[b,]) # Total signal
        CPMNSiteBasicStats[b,3] <- CPMNSiteBasicStats[b,2] / (500 + motifWidth) # Total signal per bp
        CPMNSiteBasicStats[b,4] <- sum(CPMNinsMatrix[b,(250:(250 + motifWidth))] / motifWidth) # Motif signal per bp
        CPMNSiteBasicStats[b,5] <- (sum(CPMNinsMatrix[b,200:250]) + sum(CPMNinsMatrix[b,(250 + motifWidth):(300 + motifWidth)])) / 100 # Flank signal per bp
        CPMNSiteBasicStats[b,6] <- (sum(CPMNinsMatrix[b,1:50]) + sum(CPMNinsMatrix[b,(500 + motifWidth-50):(500 + motifWidth)])) / 100 # Background signal per bp
        CPMNSiteBasicStats[b,7] <- (sum(CPMNinsMatrix[b,1:250]) + sum(CPMNinsMatrix[b,(250 + motifWidth):(500 + motifWidth)])) / 500 # Wide flank signal per bp
        CPMNSiteBasicStats[b,8] <- CPMNSiteBasicStats[b,5] / CPMNSiteBasicStats[b,6] # Flank / background
        CPMNSiteBasicStats[b,9] <- CPMNSiteBasicStats[b,4] / CPMNSiteBasicStats[b,5] # Motif / flank
        CPMNSiteBasicStats[b,10] <- CPMNSiteBasicStats[b,4] / CPMNSiteBasicStats[b,7] # Motif / wide flank
      } # end for (b in 1:numSites)
      
      #### Generate the insertion site probability vector for raw data
      CPMNInsProb <- c()
      ##
      for (c in 1:(500 + motifWidth)){
        CPMNInsProb[c] <- sum(CPMNinsMatrix[,c])
      } # end for (c in 1:(500 + motifWidth))
      ##
      CPMNTotalSignal<- sum(CPMNInsProb)
      CPMNInsProb <- CPMNInsProb / CPMNTotalSignal
      
      #### Generate null models, use BF and BH correction to parse ####
      ## Find the unique values for total signal and generate null models
      ## Generate null models based on raw signal, correct for CPM after
      uniqueTotalSignals <- unique(rawSiteBasicStats[,2])
      siteWidth <- 500 + motifWidth
      
      ## Initiate a matrix to store the mean null signal in the null model and the input signal to null model
      nullModels <- matrix(data = NA, ncol = 2, nrow = length(uniqueTotalSignals))
      colnames(nullModels) <- c("Input signal", "Avg motif signal in null model")
      CPMNnullModels <- nullModels / factorCPM
      
      ## Calculate the null models
      for (c in 1:length(uniqueTotalSignals)){
        nullVec <- generateNullFP(1000, uniqueTotalSignals[c], siteWidth, motifWidth)
        nullModels[c,1] <- uniqueTotalSignals[c]
        nullModels[c,2] <- mean(nullVec)
      } # end for (c in 1:length(uniqueTotalSignals))
      
      ## Perform a one-tailed t-test to generate a p-value for each observed motif site
      cat("Performing one-tailed t-tests on peak subset...", "\n")
      ttest <- list() # list to store the results of the t-tests
      pvalue <- c() # vector to store the p-values
      tvalue <- c() # vector to store the t-value
      
      ## Perform t-test on all sites
      cat("Performing 1-tailed ttest", "\n")
      for (d in 1:numSites){
        ## Retrieve the total signal for the current site
        currentSignal <- c(siteBasicStats[d,2])
        ## Retrieve the appropriate null model
        currentNullModel <- nullModels[which(nullModels[,1]==currentSignal),2]
        ## Perform the t-test
        ttest[[d]] <- t.test(insMatrix[d,250:(250+motifWidth)], mu=currentNullModel, alternative="less", conf.level = 0.95)
        pvalue[d] <- ttest[[d]][["p.value"]]
        tvalue[d] <- ttest[[d]][["statistic"]][["t"]]
      } # end for (d in 1:numSites)
      
      ## Get the indices of the sites that are lower than p = 0.05
      cat("Selecting p-value passing sites > 0.05", "\n")
      idxPvaluePass <- which(pvalue < 0.05)
      pvaluePass <- pvalue[idxPvaluePass]
      ppassNumSites <- length(idxPvaluePass)
      
      ## Perform bonferroni correction
      cat("Performing bonferroni correction", "\n")
      idxBFpass <- which(pvalue < (0.05/numSites))
      BFpvaluePass <- pvalue[idxBFpass]
      BFpassNumSites <- length(idxBFpass)
      
      ## Perform benjamini-hochberg correction
      BHpvalue <- p.adjust(pvalue, method = "BH")
      
      ## Data transfer to storage object and save
      parseData <- list()
      ##
      parseData$PWM <- PWM
      parseData$motifWidth <- motifWidth
      parseData$allSites <- allSites
      parseData$extSites <- extSites
      parseData$numSites <- numSites
      parseData$insMatrix <- insMatrix
      parseData$libSize <- libSize
      parseData$coverageSize <- coverageSize
      parseData$libFactor <- libFactor
      parseData$sampleTotalReads <- sampleTotalReads
      parseData$siteBasicStats <- siteBasicStats
      parseData$samplePeaks <- samplePeaks
      parseData$insStandardDeviation <- insStandardDeviation
      parseData$insMean <- insMean
      parseData$zscoreInsMatrix <- zscoreInsMatrix
      ##
      parseData$uniqueTotalSignals <- uniqueTotalSignals
      parseData$nullModels <- nullModels
      parseData$ttest <- ttest
      parseData$pvalue <- pvalue
      parseData$tvalue <- tvalue
      ##
      parseData$idxPvaluePass <- idxPvaluePass
      parseData$pvaluePass <- pvaluePass
      parseData$ppassNumSites <- ppassNumSites
      ##
      parseData$idxBFpass <- idxBFpass
      parseData$BFpvaluePass <- BFpvaluePass
      parseData$BFpassNumSites <- BFpassNumSites
      ##
      parseData$BHpvalue <- BHpvalue
      ##
      com <- paste0("footprintData$motif", a, "<- parseData")
      eval(parse(text = com))
      
    } # end for (a in 1:numMotif)
    
    ## Save the parsed footprint data
    cat("Saving finished data for", geneName, "\n")
    save(footprintData, file = dataOutPath)
    
  } # end if (numMotif == 0)
} # end if (file.exists(dataOutPath) == TRUE)

## Create the output file for snakemake
file.create(outPath)
cat("Finished parsing", "\n")



