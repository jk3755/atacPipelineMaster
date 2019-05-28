## Disable scientific notation in variables
options(scipen = 999)
## suppress warnings globally here, as they will disrupt the tryCatch block
## will need to improve this code at some point
options(warn = -1)


## Set snakemake variables
footprintDataPath <- snakemake@input[[1]]
sampleTotalReadsPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]


##
cat("Parsing and generating footprint statistics for", geneName, "from sample", sampleName, "\n")


## Set the output filepath for the Rdata object and perform a filecheck
dataOutPath <- gsub("operations/parse", "data/parsed", outPath)
dataOutPath <- gsub("parseFP.bamcopy\\d+.done", "parsedFootprintData.Rdata", dataOutPath, perl = TRUE)
cat("Output path for parsed data:", dataOutPath, "\n")


## Perform a filecheck 
if (file.exists(dataOutPath) == TRUE){
  
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
  suppressMessages(library(seqLogo))
  suppressMessages(library(ChIPpeakAnno))
  suppressMessages(library(rlist))
  
  ## Build functions
  cat("Building functions...", "\n")
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
  
  plotInsertionProb <- function(plotTitle = c(""), motifWidth, motifPWM, plotLogo = FALSE, insVector){
    
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
  
  plotInsertionHeatmap <- function(plotTitle = c(""), motifWidth, motifPWM, plotLogo = FALSE, insVector){
  #### Make heatmap for bf passing sites ####
  ## USE THIS STRUCTURE FOR HEATMAPS
  ## first, combine the signals from plus and minus strand
  heatSigs <- bfInsMatrix
  heatNumSites <- bfNumSites
  heatNumBP <- siteWidth
  heatSites <- bfSites
  
  ## scale each row individually
  for (f in 1:heatNumSites){
    maxsig <- max(heatSigs[f,])
    for (g in 1:heatNumBP){heatSigs[f,g] <- (heatSigs[f,g] / maxsig)}}
  maxsig <- 1
  ## invert signals
  for (h in 1:heatNumSites){for (i in 1:heatNumBP){heatSigs[h,i] <- (1-heatSigs[h,i])}}
  
  ## Annotate the combined sublist name which will become the tital of the heatmap plot
  heatTitle <- paste0(geneName, "_motif", a, "_numsites", heatNumSites)
  combined <- list()
  com <- paste0("combined$", heatTitle, " <- heatSigs")
  eval(parse(text = com))
  
  ##
  svgPath <- paste0(dirPath, "footprints/graphs/heatmaps/", sampleName, ".", geneName, ".", "motif", a, ".bfpeak.sites.heatmap.svg")
  svg(file = svgPath)
  cat("Saving svg footprint image at path:", svgPath, "\n")
  
  ## Margin controls
  # margin(a,b,c,d)
  # a = size of graph from top to bottom, higher value = smaller. default = 0.1
  # b = size of graph from left to right, higher value = smaller. default = 0.005
  # c = flips x axis?
  # d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
  # good settings for ATACseq = c(0.1,0.005,0.05,0.2)
  # bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values
  
  ##
  ChIPpeakAnno::featureAlignedHeatmap(combined,
                                      feature.gr = reCenterPeaks(heatSites, width = heatNumBP),
                                      upper.extreme = maxsig, # set this to control the heatmap scale
                                      annoMcols = "score",
                                      sortBy = "score",
                                      n.tile = heatNumBP,
                                      margin = c(0.1, 0.005, 0.05, 0.2),
                                      color = colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias = 0.9)(100),
                                      gp = gpar(fontsize = 10),
                                      newpage = TRUE)
  dev.off()
  } # end plotInsertionHeatmap function
  
  ###########################
  ## Make graph of the raw peak sites
  #svgPath <- paste0(dirPath, "footprints/graphs/peaks/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
  #svg(file = svgPath)
  #cat("Saving peaks footprint image at path:", svgPath, "\n")
  #plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".rawpeaks")
  #plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = insVector)
  #dev.off()
  ############################
  
  ## Get the total number of reads in the sample
  cat("Loading total sample reads from:", sampleTotalReadsPath, "\n")
  load(sampleTotalReadsPath)
  sampleTotalReads <- sampleTotalReads[6]
  cat("Found", sampleTotalReads, "total reads in current sample", "\n")
  
  ## Load the raw footprintData file
  footprintDataPath <- gsub("operations", "data", footprintDataPath)
  footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
  cat("Loading raw footprintData file from path:", footprintDataPath, "\n")
  load(footprintDataPath)
  
  ## To avoid errors, clear the list of any empty sub-lists first
  ## Remove motifs with 0 binding sites
  cat("Removing motifs that matched 0 genomic loci", "\n")
  footprintData <- list.clean(footprintData, function(footprintData) length(footprintData) == 0L, TRUE)
  ## The number of unique motifs for the current gene
  numMotif <- length(footprintData)
  cat("Found data for", numMotif, "non-0 unique motifs", "\n")
  
  
  ## If no raw data is found, skip
  if (numMotif == 0){
    cat(numMotif, "No data to analyze. Skipping", "\n")
  } else {
  
    #### Perform the parse and stats analysis 
    cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
    
    ##
    for (a in 1:numMotif){
      
      cat("Processing data for motif", a, "\n")
      
      ##
      tryCatch({
        
        ## Prepare the data from the raw footprint analysis
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
        
        
        #### Calculate basic statistics for each site ####
        siteBasicStats <- matrix(data = NA, nrow = numSites, ncol = 7)
        colnames(siteBasicStats) <- c("Site index", "Total signal", "Total signal per bp", "Motif signal per bp",
                                      "Flank signal per bp", "Background signal per bp", "Wide flank signal per bp") 
        
        
        ## Populate the basic stats matrix
        for (b in 1:numSites){
          siteBasicStats[b,1] <- b
          siteBasicStats[b,2] <- sum(insMatrix[b,])
          siteBasicStats[b,3] <- siteBasicStats[b,2] / (500 + motifWidth)
          siteBasicStats[b,4] <- sum(insMatrix[b,(250:(250+motifWidth))] / motifWidth)
          siteBasicStats[b,5] <- (sum(insMatrix[b,200:250]) + sum(insMatrix[b,(250+motifWidth):(300+motifWidth)])) / 100
          siteBasicStats[b,6] <- (sum(insMatrix[b,1:50]) + sum(insMatrix[b,(500+motifWidth-50):(500+motifWidth)])) / 100
          siteBasicStats[b,7] <- (sum(insMatrix[b,1:250]) + sum(insMatrix[b,(250+motifWidth):(500+motifWidth)])) / 500
        } # end for (b in 1:numSites)
        
        
        
        ## Normalize values in insertion matrix, z-scores
        insertionStandardDeviation <- sd(insertionMatrix)
        normalizedInsertionMatrix <- ((insertionMatrix - libraryFactor) / insertionStandardDeviation)
        
        
        ## Transfer all the data for the current motif to the storage object
        cat("Transferring all data to storage object footprintData", "\n")
        com <- paste0("footprintData$motif", idxMotif, " <- tempData")
        eval(parse(text = com))
        
        ## Update the motif index counter
        cat("Updating motif index counter", "\n")
        idxMotif <- (idxMotif + 1)
        
        
        ## Find the unique values for total signal and generate null models
        uniqueTotalSignals <- unique(siteTotalSignal)
        
        ## Initiate a matrix to store the mean null signal in the null model and the input signal to null model
        nullModels <- matrix(data = NA, ncol = 2, nrow = length(uniqueTotalSignals))
        colnames(nullModels) <- c("Input signal", "Avg motif signal in null model")
        
        ## Calculate the null models
        for (c in 1:length(uniqueTotalSignals)){
          nullVec <- generateNullFP(1000, uniqueTotalSignals[c], siteWidth, motifWidth)
          nullModels[c,1] <- uniqueTotalSignals[c]
          nullModels[c,2] <- mean(nullVec)} # end for (c in 1:length(uniqueTotalSignals))
        
        ## Perform a one-tailed t-test to generate a p-value for each observed motif site
        cat("Performing one-tailed t-tests on peak subset...", "\n")
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
        cat("Selecting p-value passing sites...", "\n")
        idxPvaluePass <- which(pvaluePeak < 0.05)
        peakPvaluePass <- pvaluePeak[idxPvaluePass]
        
        ## Perform bonferroni correction
        cat("Performing bonferroni correction...", "\n")
        idxbfPeakPass <- which(pvaluePeak < (0.05/numSites))
        bfPvaluePeakPass <- pvaluePeak[idxbfPeakPass]
        
        ## Subset the insertion matrix based on the bonferroni passing sites only
        cat("Subsetting sites based on bf corrected p-values...", "\n")
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
        
        
        # end for (a in 1:numMotif)
        # end if (numMotif == 0)
  # end if (file.exists(dataOutPath) == TRUE)
    
     
## Create the output file for snakemake
cat("Creating output file for snakemake", geneName, "\n")
file.create(outPath)





  
  
      
      ## Data transfer to storage object and save
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
      
    }, # end try
    error=function(cond){
      message(cond)
      return(NA)
    },
    finally={})
    gc()
    
  } # end for (a in 1:numMotif)
  gc()
  
  ## Finish the script and create the output file for snakemake
  save(footprintData, file = dataOutPath)
  
} # end if (file.exists(dataOutPath) == TRUE)

##
file.create(outPath)
cat("Finished!", "\n")



