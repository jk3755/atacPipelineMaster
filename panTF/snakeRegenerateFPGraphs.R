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
#suppressMessages(library(ggplot2))
#suppressMessages(library(ggpubr))
suppressMessages(library(ChIPpeakAnno))

## Set snakemake variables
cat("Setting snakemake vars...", "\n")
footprintDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]
geneName <- snakemake@wildcards[["gene"]]
dirPath <- snakemake@wildcards[["path"]]

## Load the footprintData file
cat("Loading footprintData file...", "\n")
footprintDataPath <- gsub("operations", "data", footprintDataPath)
footprintDataPath <- gsub("rawFPanalysis.bamcopy\\d+.done", "rawFootprintData.Rdata", footprintDataPath, perl = TRUE)
load(footprintDataPath)

## The number of unique motifs for the current gene
numMotif <- length(footprintData)

## Performing parsing operations
cat("Regenerating footprint graphs for", geneName, "with", numMotif, "unique motifs", "\n")
##
for (a in 1:numMotif){
  
  ## suppress warnings globally here, as they will disrupt the tryCatch block
  ## will need to improve this code at some point
  options(warn = -1)
  func <- tryCatch({
    
    com <- paste0("tempData <- footprintData$motif", a)
    eval(parse(text = com))
    motifWidth <- tempData[["motifWidth"]]
    PWM <- footprintData[["motif1"]][["PWM"]][a]
    insProfile <- tempData[["rawProfile"]]
    insVector <- insProfile[2,]
    
    ## Peaks insertion probability plot
    svgPath <- paste0(dirPath, "footprints/graphs/peaks/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
    
    if (file.exists(svgPath)){
      cat("Peaks footprint plot found at path", svgPath, "\n")
      cat("Skipping...", "\n")
      
    } else {
      
      ## Make graph of the raw peak sites
      svgPath <- paste0(dirPath, "footprints/graphs/peaks/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
      svg(file = svgPath)
      cat("Saving peaks footprint image at path:", svgPath, "\n")
      plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".rawpeaks")
      plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = insVector)
      dev.off()
      
    } # end if (file.exists(svgPath))
    
    ################################################################################
    ## bf insertion probability plot
    
    ##
    tempParse <- list()
    com <- paste0("tempParse <- footprintData$motif", a, "$parseData")
    eval(parse(text = com))
    
    ##
    bfInsMatrix <- tempParse$bfInsMatrix
    bfTotalSignal <- tempParse$bfTotalSignal
    bfProfile <- tempParse$bfProfile
    bfVector <- tempParse$bfVector
    bfSites <- tempParse$bfSites
    bfNumSites <- tempParse$bfNumSites
    
    ##
    svgPath <- paste0(dirPath, "footprints/graphs/bf/", sampleName, ".", geneName, ".", "motif", a, ".bf.sites.svg")
    ##
    if (file.exists(svgPath)){
      cat("bf corrected sites footprint plot found at", svgPath, "\n")
      cat("Skipping...", "\n")
    } else {
      svg(file = svgPath)
      cat("Saving bf footprint image at path:", svgPath, "\n")
      plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".bfsites")
      plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = bfVector)
      dev.off()
      
    } # end if (file.exists(svgPath))
    
    ################################################################################
    ## Heatmap of bf sites
    
    #### Make heatmap for bf passing sites ####
    ## USE THIS STRUCTURE FOR HEATMAPS
    ## first, combine the signals from plus and minus strand
    heatSigs <- bfInsMatrix
    heatNumSites <- bfNumSites
    heatNumBP <- siteWidth
    heatSites <- bfSites
    
    ##
    svgPath <- paste0(dirPath, "footprints/graphs/heatmaps/", sampleName, ".", geneName, ".", "motif", a, ".bfpeak.sites.heatmap.svg")
    if (file.exists(svgPath)){
      cat("heatmap plot found at", svgPath, "\n")
      cat("Skipping...", "\n")
    } else {
      
      ## scale each row individually
      for (f in 1:heatNumSites){
        maxSignal <- max(heatSigs[f,])
        for (g in 1:heatNumBP){heatSigs[f,g] <- (heatSigs[f,g] / maxSignal)}}
      maxSignal <- 1
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
                                          upper.extreme = maxSignal, # set this to control the heatmap scale
                                          annoMcols = "score",
                                          sortBy = "score",
                                          n.tile = heatNumBP,
                                          margin = c(0.1, 0.005, 0.05, 0.2),
                                          color = colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias = 0.9)(100),
                                          gp = gpar(fontsize = 10),
                                          newpage = TRUE)
      dev.off()
    } # end if (file.exists(svgPath))
    
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
  gc()   
} # end for (a in 1:numMotif)
gc()

##
file.create(outPath)
cat("Finished!", "\n")
