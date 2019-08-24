
#### TESTING ####
# bamPath <- "C:\\Users\\jsk33\\Desktop\\test\\test-REP1.u.bam"
# baiPath <- "C:\\Users\\jsk33\\Desktop\\test\\test-REP1.u.bai"
# geneName <- "MUSC"
# sampleName <- "TEST"
# sampleRep <- 1
# dirPath <- "C:\\Users\\jsk33\\Desktop\\test\\"
# snakemakeTouchPath <- "C:\\Users\\jsk33\\Desktop\\test\\touch.txt"
# gc()
#########################################################################################################

#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  
  #### Set snakemake variables ####
  bamPath <- snakemake@input[[1]]
  baiPath <- snakemake@input[[2]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  dirPath <- snakemake@wildcards[["path"]]
  snakemakeTouchPath <- snakemake@output[[1]]
  
  #### Report ####
  cat("Spooling footprint analysis", "\n")
  cat("Bam file:", bamPath, "\n")
  cat("Bai file:", baiPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Directory path:", dirPath, "\n")
  cat("Snakemake touchfile path:", snakemakeTouchPath, "\n")
  
  #### Generate functions ####
  cat("Generating functions", "\n")
  
  getAllBindingSites <- function(gene, pwmScanScore = "95%"){
    
    suppressMessages(library(MotifDb))
    cat("querying mdb", "\n")
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
      allSites <- keepStandardChromosomes(allSites, pruning.mode="coarse")
      allSites <- trim(allSites)
      allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
      for (a in 2:numUniqueMotifs){
        com <- paste0("PWM <- uniqueMotifs[[", a, "]]")
        eval(parse(text = com))
        sitesTemp <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore)
        sitesTemp <- keepStandardChromosomes(sitesTemp, pruning.mode="coarse")
        sitesTemp <- trim(sitesTemp)
        sitesTemp@elementMetadata@listData$score2 <- sitesTemp@elementMetadata@listData[["score"]] / max(sitesTemp@elementMetadata@listData[["score"]])
        allSites <- c(allSites, sitesTemp)}
      
    } else {
      
      cat("processing one unique motif", "\n")
      PWM <- uniqueMotifs[[1]]
      allSites <- Biostrings::matchPWM(PWM, genome, min.score = pwmScanScore, with.score = TRUE)
      allSites <- keepStandardChromosomes(allSites, pruning.mode="coarse")
      allSites <- trim(allSites)
      allSites@elementMetadata@listData$score2 <- allSites@elementMetadata@listData[["score"]] / max(allSites@elementMetadata@listData[["score"]])
      
    }
    cat("returning allSites", "\n")
    return(allSites)
  }
  
  generateInsMatrix <- function(bamFile, bindingSites, maxWidth, chrName, upstream = 250, downstream = 250){
    extSites <- promoters(bindingSites, upstream = upstream, downstream = (downstream + maxWidth), use.names = TRUE)
    extSites <- keepStandardChromosomes(extSites, pruning.mode="coarse")
    extSites <- trim(extSites)
    param <- ScanBamParam(which = extSites)
    bamIn <- readGAlignments(bamFile, param = param)
    grIn <- granges(bamIn)
    grIn2 <- resize(grIn, width = 1)
    grPlus <- grIn2[which(strand(grIn2) == "+")]
    grMinus <- grIn2[which(strand(grIn2) == "-")]
    grPlusShifted <- shift(grPlus, shift = 4L)
    grMinusShifted <- shift(grMinus, shift = -5L)
    shiftedInsertions <- c(grPlusShifted, grMinusShifted)
    insRLE <- coverage(shiftedInsertions)
    com <- paste0("tempRLE <- list(insRLE@listData[['", chrName, "']])")
    eval(parse(text = com))
    insRLE@listData <- tempRLE
    extSitesList <- GRangesList(extSites)
    insViews <- Views(insRLE, extSitesList)
    insMatrix <- as.matrix(insViews)
    return(insMatrix)
  }
  
  calculateBasicFootprintStatistics <- function(insMatrix, bindingSites, chrName){
    
    numSites <- length(insMatrix[,1])
    rawSiteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 18)
    
    colnames(rawSiteBasicStats) <- c("Chromosome", "Motif Start", "Motif Width",
                                     "Motif Score", "Motif Score 2",
                                     "Total signal", "Total signal per bp", "Motif signal per bp",
                                     "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp",
                                     "Flank / Background", "Motif / Flank", "Motif / Wide Flank",
                                     "Binding", "p-value", "Annotation", "Closest gene")
    
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
    
    rawSiteBasicStats[,17] <- annotations@anno@elementMetadata@listData[["annotation"]]
    rawSiteBasicStats[,18] <- annotations@anno@elementMetadata@listData[["SYMBOL"]]
    
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
        ttest <- t.test(motifSignals, mu = averageNullMotifSignal, alternative = "less", conf.level = 0.95)
        pvalue <- ttest$p.value
        rawSiteBasicStats[b,16] <- pvalue
        if (pvalue < 0.05){rawSiteBasicStats[b,15] <- 1} else {rawSiteBasicStats[b,15] <- 0}}
    }
    return(rawSiteBasicStats)
  }
  
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
  
  #### Perform a filecheck ####
  footprintDataFilepath <- paste0(dirPath, "footprints/data/", sampleName, "-REP", sampleRep, ".", geneName, ".FPdata.Rdata")
  cat("Output filepath for footprint data:", footprintDataFilepath, "\n")
  
  if (file.exists(footprintDataFilepath) == TRUE){
    cat("Footprint data file already exists, skipping operation", "\n")
  } else {
    
    #### Load libraries ####
    cat("Loading libraries", "\n")
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(Biostrings))
    suppressMessages(library(MotifDb))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(Rsamtools))
    suppressMessages(library(GenomicAlignments))
    suppressMessages(library(ChIPseeker))
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    genome <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    #### Retrieve the binding sites for the current gene ####
    cat("Retrieving binding sites", "\n")
    allSites <- getAllBindingSites(gene = geneName)
    maxWidth <- max(allSites@ranges@width)
    
    #### Set scope of analysis ####
    scope <- paste0("chr", c(1:22, "X", "Y"))
    cat("scope of analysis:", scope, "\n")
    
    #### Determine scope for current binding sites ####
    currentScope <- scope[which(scope %in% allSites@seqnames@values)]
    cat("scope of current binding sites:", currentScope, "\n")
    
    #### Generate insertion matrix for each chromosome one at a time and convert to basic footprint statistics ####
    ## Doing it this way allows the script to keep the memory usage low
    for (item in currentScope){
      cat("Generating insertion matrix for", item, "\n")
      com <- paste0("tempSites <- allSites[seqnames(allSites) == '", item, "']")
      eval(parse(text = com))
      com <- paste0("insMatrix <- generateInsMatrix(bamPath,tempSites,maxWidth,'", item, "')")
      eval(parse(text = com))
      com <- paste0(item, "SiteStatistics <- calculateBasicFootprintStatistics(insMatrix,tempSites,'", item, "')")
      eval(parse(text = com))
      rm(insMatrix)
      gc()
    }
    
    rm(allSites)
    gc()
    
    #### Merge all of the site statistics tables ####
    cat("Merging site statistics tables", "\n")
    currentTables <- as.character(paste0(currentScope, "SiteStatistics"))
    currentTablesStr <- paste(currentTables, collapse = ",")
    com <- paste0("footprintSiteStatistics <- rbind(", currentTablesStr, ")")
    eval(parse(text = com))
    
    save(footprintSiteStatistics, file = footprintDataFilepath)
    
  } # end filecheck
  
  file.create(snakemakeTouchPath)
  
}, finally = {
  file.create(snakemakeTouchPath)
})