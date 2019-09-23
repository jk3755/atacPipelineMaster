

#### Query motifDB to return binding sites of given gene ####
getAllBindingSites <- function(gene, pwmScanScore = "99%"){
  
  cat("querying motifDB", "\n")
  mdbHuman <- query(MotifDb, 'hsapiens')
  cat(length(mdbHuman), "records in database", "\n")
  geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
  cat(geneIdx, "\n")
  cat("Found", length(geneIdx), "records matching current gene", "\n")
  cat("retrieving relevant records from mdb", "\n")
  tempMotifs <- list()
  c <- 1
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1}
  cat("finding unique motifs", "\n")
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  cat("found", numUniqueMotifs, "unique motifs", "\n")
  if (numUniqueMotifs > 1){
    cat("processing more than one unique motif", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
    for (a in 2:numUniqueMotifs){
      com <- paste0("PWM <- uniqueMotifs[[", a, "]]")
      eval(parse(text = com))
      sitesTemp <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
      sitesTemp <- keepStandardChromosomes(sitesTemp, pruning.mode = "coarse")
      sitesTemp <- trim(sitesTemp)
      sitesTemp@elementMetadata@listData$score2 <- sitesTemp@elementMetadata@listData[["score"]] / max(sitesTemp@elementMetadata@listData[["score"]])
      allSites <- c(allSites, sitesTemp)}
  } else {
    cat("processing one unique motif", "\n")
    PWM <- uniqueMotifs[[1]]
    allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
    allSites <- keepStandardChromosomes(allSites, pruning.mode = "coarse")
    allSites <- trim(allSites)
    allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
  }
  
  ## Remove chrM entries, if present
  cat("Removing mitochondrial binding sites, if present", "\n")
  indexChrM <- which(seqnames(allSites) == "chrM")
  if (length(indexChrM > 0)){
    allSites <- allSites[-indexChrM]
  }
  
  cat("returning allSites", "\n")
  return(allSites)
}

#### Import a BED file as a GRanges object ####
importBED <- function(bedPath){
  bedin <- readBed(bedPath, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
  bedin <- keepStandardChromosomes(bedin, pruning.mode = "coarse")
  bedin <- trim(bedin)
  return(bedin)
}

#### Generate an insertion matrix for footprinting analysis one chromosome at a time ####
generateInsertionMatrixByChr <- function(bamFile, bindingSites, maxWidth, chrName, upstream = 250, downstream = 250){
  
  ## Expand the input binding sites to the analysis window
  extSites <- promoters(bindingSites, upstream = upstream, downstream = (downstream + maxWidth), use.names = TRUE)
  extSites <- keepStandardChromosomes(extSites, pruning.mode = "coarse")
  extSites <- trim(extSites)
  
  ## Set the parameter telling where to look in bam file
  param <- ScanBamParam(which = extSites)
  
  ## Load the relevant reads from the bam file
  bamIn <- readGAlignments(bamFile, param = param)
  
  ## Convert the GRanges and shift read ends to account for Tn5 insertion
  grIn <- granges(bamIn)
  grIn2 <- resize(grIn, width = 1)
  grPlus <- grIn2[which(strand(grIn2) == "+")]
  grMinus <- grIn2[which(strand(grIn2) == "-")]
  grPlusShifted <- shift(grPlus, shift = 4L)
  grMinusShifted <- shift(grMinus, shift = -5L)
  shiftedInsertions <- c(grPlusShifted, grMinusShifted)
  
  ## Convert to run-length encoding
  insRLE <- coverage(shiftedInsertions)
  
  ## Pull the data for the current given chromosome only
  com <- paste0("tempRLE <- list(insRLE@listData[['", chrName, "']])")
  eval(parse(text = com))
  
  insRLE@listData <- tempRLE
  extSitesList <- GRangesList(extSites)
  insViews <- Views(insRLE, extSitesList)
  
  ## Convert RLE views object to insertion matrix
  insMatrix <- as.matrix(insViews)
  
  ## Return the matrix
  return(insMatrix)
}

#### Analyze footprinting statistics given an insertion matrix ####
calculateBasicFootprintStatistics <- function(insMatrix, bindingSites, chrName){
  
  numSites <- length(insMatrix[,1])
  rawSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 19)
  
  colnames(rawSiteBasicStats) <- c("Chromosome", "Motif Start", "Motif Width",
                                   "Motif Score", "Motif Score 2",
                                   "Total signal", "Total signal per bp", "Motif signal per bp",
                                   "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                   "Flank / Background", "Motif / Flank", "Motif / Wide Flank",
                                   "Binding", "p-value", "z-score", "Annotation", "Closest gene")
  
  rawSiteBasicStats[,1] <- chrName
  rawSiteBasicStats[,2] <- as.numeric(bindingSites@ranges@start)
  rawSiteBasicStats[,3] <- as.numeric(bindingSites@ranges@width)
  rawSiteBasicStats[,4] <- as.numeric(bindingSites@elementMetadata@listData[["score"]])
  rawSiteBasicStats[,5] <- as.numeric(bindingSites@elementMetadata@listData[["score2"]])
  
  ## Annotate
  annotations <- annotatePeak(
    peak = bindingSites,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    level = "gene",
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
    annoDb = "org.Hs.eg.db",
    addFlankGeneInfo = FALSE,
    flankDistance = 5000,
    sameStrand = FALSE,
    ignoreOverlap = FALSE,
    ignoreUpstream = FALSE,
    ignoreDownstream = FALSE,
    overlap = "TSS",
    verbose = TRUE)
  
  rawSiteBasicStats[,18] <- annotations@anno@elementMetadata@listData[["annotation"]]
  rawSiteBasicStats[,19] <- annotations@anno@elementMetadata@listData[["SYMBOL"]]
  
  for (b in 1:numSites){
    tempWidth <- as.numeric(rawSiteBasicStats[b,3])
    tempTotalWidth <- as.numeric(tempWidth + 500)
    tempTotalSignal <- as.numeric(sum(insMatrix[b,1:tempTotalWidth]))
    tempTotalSignalPerBP <- as.numeric(tempTotalSignal/tempTotalWidth)
    tempMotifSignalPerBP <- as.numeric(sum(insMatrix[b,(250:(250 + tempWidth))] / tempWidth))
    tempFlankSignalPerBP <- as.numeric((sum(insMatrix[b,200:250])+sum(insMatrix[b,(250+tempWidth):(300+tempWidth)])) / 100)
    tempBackgroundSignalPerBP <- as.numeric((sum(insMatrix[b,1:50])+sum(insMatrix[b,(500+tempWidth-50):(500+tempWidth)])) / 100)
    tempWideFlankSignalPerBP <- as.numeric((sum(insMatrix[b,1:250])+sum(insMatrix[b,(250+tempWidth):(500+tempWidth)])) / 500)
    ##
    rawSiteBasicStats[b,6] <- as.numeric(tempTotalSignal)
    rawSiteBasicStats[b,7] <- as.numeric(tempTotalSignalPerBP)
    rawSiteBasicStats[b,8] <- as.numeric(tempMotifSignalPerBP)
    rawSiteBasicStats[b,9] <- as.numeric(tempFlankSignalPerBP)
    rawSiteBasicStats[b,10] <- as.numeric(tempBackgroundSignalPerBP)
    rawSiteBasicStats[b,11] <- as.numeric(tempWideFlankSignalPerBP)
    rawSiteBasicStats[b,12] <- as.numeric(tempFlankSignalPerBP / tempBackgroundSignalPerBP)
    rawSiteBasicStats[b,13] <- as.numeric(tempMotifSignalPerBP / tempFlankSignalPerBP)
    rawSiteBasicStats[b,14] <- as.numeric(tempMotifSignalPerBP / tempWideFlankSignalPerBP)
    
    if (tempTotalSignal == 0){
      
    } else {
      averageNullMotifSignal <- generateNullFP(1000, tempTotalSignal, tempTotalWidth, tempWidth)
      motifSignals <- c(insMatrix[b,(250:(250 + tempWidth))])
      
      result = tryCatch({
        ttest <- t.test(motifSignals, mu = averageNullMotifSignal, alternative = "less", conf.level = 0.95)
        pvalue <- ttest$p.value
        rawSiteBasicStats[b,16] <- pvalue
        rawSiteBasicStats[b,17] <- qnorm(pvalue)
        if (pvalue < 0.05){rawSiteBasicStats[b,15] <- 1} else {rawSiteBasicStats[b,15] <- 0}
      }, warning = function(w) {
      }, error = function(e) {
      }, finally = {
      })}
  }
  
  #### Convert to a dataframe before return ####
  ## This allows you to store different data types
  cat("Transferring data to dataframe", "\n")
  rawFootprintStats <- data.frame(
    chr = as.character(rawSiteBasicStats[,1]),
    start = as.numeric(rawSiteBasicStats[,2]),
    width = as.numeric(rawSiteBasicStats[,3]),
    score = as.numeric(rawSiteBasicStats[,4]),
    score2 = as.numeric(rawSiteBasicStats[,5]),
    totalSignal = as.numeric(rawSiteBasicStats[,6]),
    totalSignalPerBP = as.numeric(rawSiteBasicStats[,7]),
    motifSignalPerBP = as.numeric(rawSiteBasicStats[,8]),
    flankSignalPerBP = as.numeric(rawSiteBasicStats[,9]),
    backgroundSignalPerBP = as.numeric(rawSiteBasicStats[,10]),
    wideFlankSignalPerBP = as.numeric(rawSiteBasicStats[,11]),
    flankAccessibility = as.numeric(rawSiteBasicStats[,12]),
    footprintDepth = as.numeric(rawSiteBasicStats[,13]),
    wideAccessibility = as.numeric(rawSiteBasicStats[,14]),
    binding = as.numeric(rawSiteBasicStats[,15]),
    pvalue = as.numeric(rawSiteBasicStats[,16]),
    zscore = as.numeric(rawSiteBasicStats[,17]),
    annotation = as.character(rawSiteBasicStats[,18]),
    symbol = as.character(rawSiteBasicStats[,19])
  )
  cat("Returning dataframe", "\n")
  return(rawFootprintStats)
}

#### Generate a null model of a TF binding footprint for hypothesis testing ####
generateNullFP <- function(iterations, inputSignal, analysisWidth, motifWidth){
  averages <- c()
  for (a in 1:iterations){
    null <- c(1:(analysisWidth))
    null <- c(as.vector(rmultinom(1, size=inputSignal, prob=rep(1, length(null)))))
    motifStart <- ((analysisWidth - motifWidth)/2)
    motifEnd <- (motifStart + motifWidth)
    motifAvg <- (sum(null[motifStart:motifEnd])) / motifWidth
    averages[a] <- motifAvg
  }
  return(mean(averages))
}

#### Merge two or more GRanges objects ####
getGRangesUnion <- function(grList = list()){
  numItems <- length(grList)
  cat("Processing the union of", numItems, "GRanges objects", "\n")
  ## Unlist the GRanges
  for (a in 1:numItems){
    com <- paste0("tempGR", a, " <- grList[[", a, "]]")
    eval(parse(text = com))}
  ## Merge the first two
  cat("Merging first two GRanges", "\n")
  unionGR <- union(tempGR1, tempGR2)
  ## Merge the rest
  cat("Merging remaining GRanges", "\n")
  if (numItems > 2){
    for (b in 3:numItems){
      com <- paste0("unionGR <- GenomicRanges::union(unionGR, tempGR", b, ")")
      eval(parse(text = com))}}
  ## Return
  return(unionGR)}

#### Write a GRanges object to a 3-column BED file ####
writeGRangesToBED <- function(inputGR, filePath){
  
  ## BED files are 0-coordinate based, so subtract start by 1
  cat("Converting GRanges object to BED format", "\n")
  dataframeGR <- data.frame(
    chrom = seqnames(inputGR),
    chromStart = start(inputGR)-1,
    chromEnd = end(inputGR)
  )
  ## Write the output file
  cat("Writing BED output file at path:", filePath, "\n")
  write.table(
    dataframeGR,
    file = filePath,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
}

#### Write a GRanges object to a 3-column BED file ####
writeGRangesToBEDwithStrand <- function(inputGR, filePath){
  
  ##
  numRecords <- length(inputGR)
  
  ## BED files are 0-coordinate based, so subtract start by 1
  ## Must keep name and score columns for readBed to work with strand
  cat("Converting GRanges object to BED format", "\n")
  dataframeGR <- data.frame(
    chrom = seqnames(inputGR),
    chromStart = start(inputGR)-1,
    chromEnd = end(inputGR),
    name = c(1:numRecords),
    score = rep(1, time = numRecords),
    strand = as.character(strand(inputGR))
  )
  ## Write the output file
  cat("Writing BED output file at path:", filePath, "\n")
  write.table(
    dataframeGR,
    file = filePath,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = F
  )
}

#### Generate motif insertion probability graph ####
plotInsertionSiteProbability <- function(insertionMatrix, motifStart, motifWidth, plotTitle, svgOutPath){
  
  ## Get matrix data
  matrixTotalSignal <- sum(insertionMatrix)
  columnTotalSignal <- colSums(insertionMatrix)
  rowTotalSignal <- rowSums(insertionMatrix)
  insertionProbability <- (columnTotalSignal / matrixTotalSignal) * 100
  maxInsertionProbability <- max(insertionProbability)
  totalBP <- length(columnTotalSignal)
  flankBP <- ((totalBP - motifWidth) / 2 )
  
  ## Set the plot labels
  xLabel = "Distance to motif (bp)"
  yLabel = "Tn5 fragmentation probability"
  
  ## Set the plot axis limits
  xLimit <- c(0, totalBP + 1)
  yLimit <- c(0, maxInsertionProbability * 1.5)
  
  #### Set up the plot viewports ####
  grid.newpage()
  
  ## Set the outer bounds of the figure
  vpRegion <- plotViewport(margins=c(5.1, 5.1, 4.1, 2.1), name="plotRegion")
  pushViewport(vpRegion)
  
  ## Setup the graph plotting area
  vpGraph <- viewport(
    y = .4, 
    height = .8,
    xscale = xLimit,
    yscale = yLimit,
    name="graphRegion")
  
  pushViewport(vpGraph)
  
  ## Add the insertion probability data line
  # lwd = line width
  # col = line color
  grid.lines(
    x = 1:totalBP,
    y = insertionProbability,
    default.units = "native",
    gp = gpar(lwd = 1, col = "darkred"))
  
  ## Add the x-axis and tick marks
  # at = is a numeric vector with the x-axis locations for tick marks
  # label = specify labels on tick marks
  grid.xaxis(at = 
               c(seq(1, flankBP, length.out = 3),
                 flankBP + seq(1, motifWidth),
                 flankBP + motifWidth + seq(1, flankBP, length.out = 3)),
             label = c(-(flankBP + 1 - seq(1, flankBP + 1, length.out = 3)),
                       rep("", motifWidth),
                       seq(0, flankBP, len = 3)))
  
  ## Add the y-axis
  grid.yaxis()
  
  ## Adds the dashed lines bounding the motif
  grid.lines(
    x = c(flankBP, flankBP, 0),
    y = c(0, max(insertionProbability), yLimit[2]),
    default.units = "native",
    gp = gpar(lty = 2))
  grid.lines(
    x = c(flankBP + motifWidth + 1, flankBP + motifWidth + 1, totalBP),
    y = c(0, maxInsertionProbability, yLimit[2]),
    default.units = "native",
    gp = gpar(lty = 2))
  
  #### Add the plot title and x and y axis labels
  upViewport()
  vpTitle <- viewport(
    y = .9,
    height = .2,
    xscale = c(0, totalBP + 1),
    name = "title")
  
  pushViewport(vpTitle)
  upViewport()
  
  vpLegend <- viewport(
    x = 0.5,
    y = 0.5,
    width = convertX(unit(1, "lines"), unitTo = "npc"),
    height = convertY(unit(1, "lines"), unitTo = "npc"),
    just = c("right", "top"), name = "legendWraper")
  pushViewport(vpLegend)
  upViewport()
  
  grid.text(
    plotTitle,
    y = unit(1, "npc")-convertY(unit(1, "lines"), unitTo = "npc"),
    gp = gpar(cex = 1.2, fontface = "bold"))
  upViewport()
  
  grid.text(xLabel, y = unit(1, 'lines'))
  grid.text(yLabel, x = unit(1, 'line'), rot = 90)
  
  #### Save the image as an svg
  grid.export(svgOutPath)
  dev.off()
}

#### Generate a motif-aligned signal heatmap ####
plotMotifAlignedHeatmap <- function(insertionMatrix, bindingSites, plotTitle, svgOutPath){
  
  #### Margin controls ####
  # margin(a,b,c,d)
  # a = size of graph from top to bottom, higher value = smaller. default = 0.1
  # b = size of graph from left to right, higher value = smaller. default = 0.005
  # c = flips x axis?
  # d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
  # good settings for ATACseq = c(0.1,0.005,0.05,0.2)
  # bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values
  # upper.extreme controls to heatmap scale
  
  ## Get the matrix data
  numSites <- length(insertionMatrix[,1])
  numBP <- length(insertionMatrix[1,])
  
  ## Create list structure required to run the featureAlignedHeatmap function
  plotData <- list()
  com <- paste0("plotData$", plotTitle, " <- insertionMatrix")
  eval(parse(text = com))
  
  ## Make the plot
  ChIPpeakAnno::featureAlignedHeatmap(
    plotData,
    feature.gr = reCenterPeaks(bindingSites, width = numBP),
    upper.extreme = 1,
    annoMcols = "score",
    sortBy = "score",
    n.tile = numBP,
    margin = c(0.1, 0.005, 0.05, 0.2),
    color = colorRampPalette(c("white", "red"), bias = 0.5)(100),
    gp = gpar(fontsize = 10),
    newpage = TRUE)
  grid.export(svgOutPath)
  dev.off()
}

#### Generate a seqbias correction model ####
generateTn5BiasCorrectionModel <- function(refFastaPath, bamPath, refBEDpath, biasedPlotPath, correctedPlotPath, biasModelPath){
  
  # refFastaPath - path to the reference fasta file
  # bamPath - path to the input bam file
  # refGR - a GRanges object containing the intervals to examine for the model
  # biasedPlotPath - output path for the biased insertion freqs
  # correctedPlotPath - output path for the corrected insertion freqs
  # biasModelPath - output path to save the correction model
  
  #### Load the BED interval ####
  cat("Loading interval from BED file", "\n")
  refGR <- importBED(refBEDpath)
  # need to drop * strand level
  refGR@strand@values <- droplevels(refGR@strand@values)
  
  #### Load the reference fasta file ####
  cat("Loading reference fasta", "\n")
  refFasta <- FaFile(refFastaPath)
  open.FaFile(refFasta)
  refSeqIn <- scanFa(refFasta, refGR)

  #### Get the reverse stand indices and reverse complement ####
  cat("Getting negative strand indices", "\n")
  negIdx <- as.logical(refGR@strand == '-')
  refSeqIn[negIdx] <- reverseComplement(refSeqIn[negIdx])
  
  #### Count the reads in the given interval ####
  cat("Counting reads in given intervals", "\n")
  readCounts <- count.reads(bamPath, refGR, binary = T)
  readFreqs <- kmer.freq(refSeqIn, readCounts)
  
  #### Plot the biased sequence plot ####
  cat("Generating biased plot", "\n")
  svg(biasedPlotPath)
  P <- qplot(
            x = pos,
            y = freq,
            ylim = c(0, 0.5),
            color = seq,
            data = readFreqs,
            geom = "line"
            )
  P <- P + facet_grid( seq ~ . )
  print(P)
  dev.off()
  
  #### Train the compensation model ####
  cat("Training compensation model", "\n")
  seqBiasModel <- seqbias.fit(
                              refFastaPath,
                              bamPath,
                              L = 20,
                              R = 20
                              )
  
  #### Predict sequence bias ####
  cat("Predicting insertion bias", "\n")
  ## The seqbias prediction must be given the underlying sequences you want to correct
  ## The model can be trained on chr1, etc.
  seqBiasPrediction <- seqbias.predict(
                                      seqBiasModel,
                                      refGR
                                      )
  
  #### Adjust sequence bias ####
  cat("Adjusting insertion counts", "\n")
  adjustedCounts <- mapply(
                          FUN = `/`,
                          readCounts,
                          seqBiasPrediction,
                          SIMPLIFY = F
                          )
  
  #### Check and plot adjustment ####
  cat("Adjusting insertion frequencies", "\n")
  adjustedReadFreqs <- kmer.freq(
                                refSeqIn,
                                adjustedCounts
                                )
  
  #### Generate and save the corrected freq plot ####
  cat("Generating corrected plot", "\n")
  svg(correctedPlotPath)
  P <- qplot( 
            x = pos,
            y = freq,
            ylim = c(0.0, 0.5),
            color = seq,
            data = adjustedReadFreqs,
            geom = "line"
            )
  P <- P + facet_grid( seq ~ . )
  print(P)
  dev.off()
  
  #### Return and save the model ####
  cat("Saving seqbias model", "\n")
  seqbias.save(seqBiasModel, biasModelPath)
  #return(seqBiasModel)
}

#### Generate bias uncorrected null footprint model ####
generateFootprintNullModelUncorrected <- function(inputSiteGR, iterations = 10000, inputSignal){
  
  #### Get the width of the current motif ####
  motifWidth <- inputSiteGR@ranges@width
  analysisWidth <- motifWidth + 500
  
  #### Generate the random insertion matrix with no correction ####
  nullMatrixUncorrected <- matrix(data = NA, nrow = iterations, ncol = analysisWidth)
  
  for (a in 1:iterations){
    nullMatrixUncorrected[a,] <- rmultinom(
                                          n = 1,
                                          size = inputSignal,
                                          prob = rep(1, analysisWidth))}
  
  #### Collect the average signals ####
  flank1AvgSignalUncorrected <- c()
  flank2AvgSignalUncorrected <- c()
  motifAvgSignalUncorrected <- c()
  
  for (b in 1:iterations){
    
    flank1AvgSignalUncorrected[b] <- mean(nullMatrixUncorrected[b , 230:249 ])
    flank2AvgSignalUncorrected[b] <- mean(nullMatrixUncorrected[b , (251 + motifWidth):(270 + motifWidth) ])
    motifAvgSignalUncorrected[b]  <- mean(nullMatrixUncorrected[b , 250:(250 + motifWidth) ])
    
  }
  
  #### Transfer data to df for return ####
  nullDataUncorrected <- data.frame(
                                    flank1AvgSignal = flank1AvgSignalUncorrected,
                                    flank2AvgSignal = flank2AvgSignalUncorrected,
                                    motifAvgSignal = motifAvgSignalUncorrected,
                                    flankAvgSignal = ((flank1AvgSignalUncorrected + flank2AvgSignalUncorrected) / 2),
                                    signalDiff = ((flank1AvgSignalUncorrected + flank2AvgSignalUncorrected) / 2) - motifAvgSignalUncorrected
                                    )
  
  return(nullDataUncorrected)
}

#### Generate Tn5 bias-corrected null footprint model ####
generateFootprintNullModelBiasCorrected <- function(seqbiasModel, inputSiteGR, iterations = 10000, inputSignal){
  
  #### Get the width of the current motif ####
  motifWidth <- inputSiteGR@ranges@width
  analysisWidth <- motifWidth + 500
  
  #### Expand the input binding site GRanges to the analysis window ####
  extSites <- promoters(inputSiteGR, upstream = 250, downstream = (250 + motifWidth), use.names = TRUE)
  extSites <- keepStandardChromosomes(extSites, pruning.mode = "coarse")
  extSites <- trim(extSites)
  
  #### Pull the underlying sequence and add to the GRanges ####
  extSites@elementMetadata@listData[["string"]] <- getSeq(genome, extSites)
  
  #### Generate the prediction model ####
  seqBiasPrediction <- seqbias.predict(seqbiasModel, extSites)
  ## Invert the prediction to be able to use as a probability vector
  inverseSeqBiasPrediction <- 1 / seqBiasPrediction[[1]]

  #### Generate the random insertion matrix with correction ####
  nullMatrixCorrected <- matrix(data = NA, nrow = iterations, ncol = analysisWidth)
  
  for (a in 1:iterations){
    nullMatrixCorrected[a,] <- rmultinom(
      n = 1,
      size = inputSignal,
      prob = inverseSeqBiasPrediction)
  }
  
  #### Collect the average signals ####
  flank1AvgSignalCorrected <- c()
  flank2AvgSignalCorrected <- c()
  motifAvgSignalCorrected <- c()
  
  for (b in 1:iterations){
    
    flank1AvgSignalCorrected[b] <- mean(nullMatrixCorrected[b , 230:249 ])
    flank2AvgSignalCorrected[b] <- mean(nullMatrixCorrected[b , (251 + motifWidth):(270 + motifWidth) ])
    motifAvgSignalCorrected[b]  <- mean(nullMatrixCorrected[b , 250:(250 + motifWidth) ])
    
  }
  
  #### Transfer data to df for return ####
  nullDataCorrected <- data.frame(
    flank1AvgSignal = flank1AvgSignalCorrected,
    flank2AvgSignal = flank2AvgSignalCorrected,
    motifAvgSignal = motifAvgSignalCorrected,
    flankAvgSignal = ((flank1AvgSignalCorrected + flank2AvgSignalCorrected) / 2),
    signalDiff = ((flank1AvgSignalCorrected + flank2AvgSignalCorrected) / 2) - motifAvgSignalCorrected
  )
  return(nullDataCorrected)
}
