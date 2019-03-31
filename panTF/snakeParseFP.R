
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
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ChIPpeakAnno))

## Set snakemake variables
cat("Setting snakemake vars...", "\n")
footprintDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

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

## Load the footprintData file
cat("Loading footprintData file...", "\n")
footprintDataPath <- gsub("operations", "data/raw", footprintDataPath)
footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
load(footprintDataPath)

## The number of unique motifs for the current gene
numMotif <- length(footprintData)

## Performing parsing operations
cat("Parsing footprint data for gene", geneName, "with", numMotif, "unique motifs", "\n")
##
for (a in 1:numMotif){
    
    ## Prepare the data
    com <- paste0("tempData <- footprintData$motif", a)
    eval(parse(text = com))
    peakSites <- tempData[["sites"]]
    numSites <- length(tempData[["insMatrix"]][,1])
    siteWidth <- length(tempData[["insMatrix"]][1,])
    motifWidth <- tempData[["motifWidth"]]
    PWM <- footprintData[["motif1"]][["PWM"]][a]
    insMatrix <- tempData[["insMatrix"]]
    insProfile <- tempData[["rawProfile"]]
    insVector <- insProfile[2,]
    siteTotalSignal <- c()
    
    ## Make graph of the raw peak sites
    svgPath <- paste0(dirpath, "footprints/graphs/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
    svg(file = svgPath)
    cat("Saving peaks footprint image at path:", svgPath, "\n")
    plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".rawpeaks")
    plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = insVector)
    dev.off()
    
    ## Calculate total signal for each site
    for (b in 1:numSites){
      siteTotalSignal[b] <- sum(tempData[["insMatrix"]][b,])} # end for (b in 1:numSites)
    
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
    
    ## Make a plot for the bf passing sites
    svgPath <- paste0(dirpath, "footprints/graphs/", sampleName, ".", geneName, ".", "motif", a, ".bf.sites.svg")
    svg(file = svgPath)
    cat("Saving peaks footprint image at path:", svgPath, "\n")
    plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".bfsites")
    plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = bfVector)
    dev.off()
      
    ## suppress warnings globally here, as they will disrupt the tryCatch block
    ## will need to improve this code at some point
    options(warn = -1)
      
    func <- tryCatch({
        
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
    svgPath <- paste0(dirpath, "footprints/heatmaps/", samplename, ".", genename, ".", "motif", x, ".bfpeak.sites.heatmap.svg")
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

    
    
    
        
        ## Data transfer to storage object and save
        parsedSitesInfo <- list()
        parsedSitesInfo$peaks <- gr_narrowPeak
        parsedSitesInfo$motif <- x
        parsedSitesInfo$PWM <- PWM
        parsedSitesInfo$motifWidth <- motifWidth
        parsedSitesInfo$libraryFactor <- libraryFactor
        ## Whole genome sites
        parsedSitesInfo$genomeSignals <- genomeSignals
        parsedSitesInfo$genomeSites <- genomeSites
        parsedSitesInfo$numGenomeSites <- numGenomeSites
        parsedSitesInfo$combinedGenomeSignal <- combinedGenomeSignal
        parsedSitesInfo$genomeSignalTotals <- genomeSignalTotals
        parsedSitesInfo$motifGenomeSignalTotals <- motifGenomeSignalTotals
        ## Raw peak overlapping sites
        parsedSitesInfo$peakSignals <- peakSignals
        parsedSitesInfo$peakSites <- peakSites
        parsedSitesInfo$numPeakSites <- numPeakSites
        parsedSitesInfo$combinedPeakSignal <- combinedPeakSignal
        parsedSitesInfo$peakProfile <- peakProfile
        parsedSitesInfo$peakSignalTotals <- peakSignalTotals
        parsedSitesInfo$motifPeakSignalTotals <- motifPeakSignalTotals
        ## BF corrected p-value passing peak sites
        parsedSitesInfo$bfPassPeakSignals <-bfPassPeakSignals
        parsedSitesInfo$bfPassPeakSites <- bfPassPeakSites
        parsedSitesInfo$numbfPassPeakSites <- numbfPassPeakSites
        parsedSitesInfo$combinedbfPassPeakSignal <- combinedbfPassPeakSignal
        parsedSitesInfo$bfPassPeakProfile <- bfPeakProfile
        parsedSitesInfo$bfPassPeakSignalTotals <- bfPassPeakSignalTotals
        parsedSitesInfo$bfPassPeakMotifPeakSignalTotals <- bfPassPeakMotifSignalTotals
        ## bf sites heatmap info
        parsedSitesInfo$combined <- combined
        parsedSitesInfo$heatsites <- heatsites
        parsedSitesInfo$heatnumbp <- heatnumbp
        ##
        save(parsedSitesInfo, file = info_path)
        
      }, # end try
      error=function(cond){
        message(cond)
        return(NA)
      },
      finally={})
    
} # end for (a in 1:numMotif)


cat("Finished...", "\n")
file.create(outpathdone)