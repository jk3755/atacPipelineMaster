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
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]
mergedpath <- snakemake@input[[4]]
peakpath <- snakemake@input[[5]]
outpathdone <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
genename <- snakemake@wildcards[["gene"]]
dirpath <- snakemake@wildcards[["path"]]

##
cat("Importing peaks file...", "\n")
file_narrowPeak = peakpath
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
gr_narrowPeak <- import(file_narrowPeak, format = "BED",
                        extraCols = extraCols_narrowPeak)

##
cat("Building functions", "\n")
#
calcLibFactor <- function(mt, bampath, baipath){
  
  mt <- keepStandardChromosomes(mt, pruning.mode="coarse")
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
  return(libFactor)
} # end calcLibFactor
#
buildProfile <- function(sigs, libFactor, upstream, wid, downstream){
  
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
  profileInfo <- list()
  profileInfo$segmentation <- Profile.seg
  profileInfo$profile <- Profile
  return(profileInfo)
} # end buildProfile
#
pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
} # end pwm2pfm
#
generateNullFP <- function(n, t, s, m){
  #This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
  #- Consider the total signal (number of insertions) at each specific ~200 bp locus
  #- Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
  #- Use published or experimentally derived models of Tn5 sequence specific insertion bias
  #- For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
  #- Generate the null model graph by weighted random residstribution of the total observed signal at that site
  #- Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
  # These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
  # n = iterations
  # t = unique values for total signal
  # s = total bp in region of interest (flank + motif)
  # m = motif width
  
  # declare vector of size n to store average motif signal values
  averages <- c()
  
  # generate the null models and calculate motif averages
  for (a in 1:n){
    
    # declare the null vector
    null <- c(1:(s))
    
    # randomly distribute the total signal
    # size = the number of values to distribute
    # prob = probability of each site
    # length = length of the generated vector
    null <- c(as.vector(rmultinom(1, size=t, prob=rep(1, length(null)))))
    
    # calculate the average signal in the motif region
    b <- (s - 100 - m) # start of motif
    c <- (s - 100) # end of motif
    avg <- sum(null[b:c])
    avg <- (avg/m)
    
    # store the values
    averages[a] <- avg
    
    ### plotting
    # make the density plot
    # generate the density distribution
    # adjust can be modified to adjust the kernel estimation bandwidth
    #d <- density(averages, adjust=1, from=0)
    # make the plot
    # main can be specified to make a chart title
    # make the title string
    #title <- paste0("Motif signal means", "\n", "n = ", n, ", input signal = ", t, ", motif size = ", m)
    #plot(d, main=title)
    ## check for uniform distribution with QQ plot
    #qqplot(averages, runif(1000))
    #abline(0,1)
    
  } # end for (a in 1:n)
  return(averages)
} # end generateNullFP function


##
cat("Parsing binding sites for gene", genename, " from sample", samplename, "\n")
load(sitespath)
num_motifs <- length(bindingSites)
cat("Found ", num_motifs, " motifs", "\n")


for (x in 1:num_motifs){
  
  info_path <- paste0(dirpath, "footprints/parsed/", samplename, ".", genename, ".", "motif", x, ".info.Rdata")
  cat("Output path for parsed binding sites object: ", info_path, "\n")
  
  if (file.exists(info_path) == TRUE){
    
    cat("File already exists, skipping...", "\n")
    next
    
  } else {
    
    cat("Parsed file not found, generating", "\n")
    ##
    cat("Setting parameters...", "\n")
    motif_score <- "95%"
    upstream <- 100
    downstream <- 100
    genome <- Hsapiens
    anchor <- c("cut site")
    
    ##
    cat("Loading data...", "\n")
    mergein <- gsub("merged.done.txt", paste0("motif", x, ".merged.Rdata"), mergedpath)
    mergein <- gsub("operations", "data/merged", mergein)
    cat("Loading merged signal file from path:", mergein, "\n")
    load(mergein)
    genomeSignals <- merged_signals[["signal"]]
    bpTotal <- length(genomeSignals[["+"]][1,])
    genomeSites <- merged_signals$bindingSites
    # create an index for the genome-wide sites
    genomeSites@elementMetadata@listData$index <- 1:length(genomeSites@ranges@start)
    PWM <- bindingSites[[x]][["PWM"]]
    motifWidth <- length(PWM[1,])
    scope <- seqlevels(genomeSites)
    
    ## Calculate site signal totals for genome sites
    # for each site, calculate the total signal and motif signal
    cat("Calculating signal at each site...", "\n")
    numGenomeSites <- length(genomeSignals[["+"]][,1])
    genomeSignalTotals <- c()
    motifGenomeSignalTotals <- c()
    for (i in 1:numGenomeSites){
      genomeSignalTotals[i] <- sum(genomeSignals[["+"]][i,] + genomeSignals[["-"]][i,])
      motifGenomeSignalTotals[i] <- sum(genomeSignals[["+"]][i,100:(100+motifWidth)] + genomeSignals[["-"]][i,100:(100+motifWidth)])}
    
    ## Generate combined signal, keep this list structure, as it is necessary for downstream featureAlignedHeatmap generation
    combinedGenomeSignal <- list()
    combinedGenomeSignal$signal <- matrix(data = NA, nrow = numGenomeSites, ncol = bpTotal)
    for (i in 1:numGenomeSites){combinedGenomeSignal$signal[i,] <- genomeSignals[["+"]][i,] + genomeSignals[["-"]][i,]}
    
    ## Calculate libFactor
    cat("Calculating libFactor...", "\n")
    libraryFactor <- calcLibFactor(genomeSites, bampath, baipath)
    
    ## First parse - binding sites in peaks only
    peakSites <- subsetByOverlaps(genomeSites, gr_narrowPeak)
    idx <- peakSites@elementMetadata@listData[["index"]]
    peakSignals <- list()
    peakSignals$"+" <- genomeSignals[["+"]][idx,]
    peakSignals$"-" <- genomeSignals[["-"]][idx,]
    peakSites <- genomeSites[idx]
    peakSignalTotals <- genomeSignalTotals[idx]
    motifPeakSignalTotals <- motifGenomeSignalTotals[idx]
    numPeakSites <- length(peakSignals[["+"]][,1])
    
    ## Generate combined signals for peak sites
    combinedPeakSignal <- list()
    combinedPeakSignal$signal <- matrix(data = NA, nrow = numPeakSites, ncol = bpTotal)
    for (i in 1:numPeakSites){combinedPeakSignal$signal[i,] <- peakSignals[["+"]][i,] + peakSignals[["-"]][i,]}
    
    ## Build profile for peak sites
    cat("Building profile...", "\n")
    peakProfile <- buildProfile(peakSignals, libraryFactor, upstream, motifWidth, downstream)
    prof <- peakProfile$profile
    ## Make graph
    svg_path <- paste0(dirpath, "footprints/graphs/", samplename, ".", genename, ".", "motif", x, ".peak.sites.svg")
    svg(file = svg_path) # set the filepath for saving the svg figure
    cat("Saving peaks footprint image at path:", svg_path, "\n")
    PWMin <- pwm2pfm(PWM)
    cat("Plotting graph...", "\n")
    ATACseqQC:::plotFootprints(prof,
                               Mlen=motifWidth, motif=PWMin)
    dev.off()
    
    ## Implement null model to select peak sites
    cat("Generating null models for raw peak sites...", "\n")
    uniqueSignals <- unique(peakSignalTotals)
    nullSims <- list()
    numUnique <- length(uniqueSignals)
    nullMean <- c()
    ##
    for (i in 1:numUnique){
      n <- 1000
      t <- uniqueSignals[i]
      s <- bpTotal
      m <- motifWidth
      nullSims[[i]] <- generateNullFP(n,t,s,m)
      nullMean[i] <- mean(nullSims[[i]])}
    
    
    ## Perform a one-tailed t-test to generate a p-value for each observed motif site
    cat("Performing one-tailed t-tests on peak subset...", "\n")
    ttestPeak <- list() # list to store the results of the t-tests
    pvaluePeak <- c() # vector to store the p-values
    tvaluePeak <- c() # vector to store the t-value
    ## Need to combine the signals into one matrix (+/- strand)
    ## use this list structure as it is necessary for downstream featureAlignedHeatmaps
    #combinedGenomeSignal <- list()
    #combinedGenomeSignal$signal <- matrix(data = NA, nrow = numPeakSites, ncol = bpTotal)
    # combine the signals
    #for (i in 1:numPeakSites){combinedGenomeSignal$signal[i,] <- peakSignals[["+"]][i,] + peakSignals[["-"]][i,]}
    
    ## 
    for (i in 1:numPeakSites){
      currentSignal <- c(combinedPeakSignal$signal[i,])
      # retrieve the appropriate null model
      currentNullModel <- nullSims[[which(uniqueSignals==sum(currentSignal))]]
      # do the t-test
      ttestPeak[[i]] <- t.test(currentSignal[100:(100+motifWidth)], mu=mean(currentNullModel), alternative="less", conf.level = 0.95)
      pvaluePeak[i] <- ttestPeak[[i]][["p.value"]]
      tvaluePeak[i] <- ttestPeak[[i]][["statistic"]][["t"]]
    }
    
    ## get the indices of the sites that are lower than p = 0.05
    cat("Selecting p-value passing sites...", "\n")
    idxPeakPvaluePass <- which(pvaluePeak < 0.05)
    peakPvaluePass <- pvaluePeak[idxPeakPvaluePass]
    
    ## bonferroni correction
    cat("Performing bonferroni correction...", "\n")
    idxbfPeakPass <- which(pvaluePeak < (0.05/numPeakSites))
    bfPeakPass <- pvaluePeak[idxbfPeakPass]
    
    ## Subset the peak sites based on bf-corrected p-values
    cat("Subsetting binding sites based on bf corrected p-values...", "\n")
    
    
    ## suppress warnings globally here, as they will disrupt the tryCatch block
    ## will need to improve this code at some point
    options(warn=-1)
    
    func <- tryCatch({
    
      bfPassPeakSignals <- list()
      bfPassPeakSignals$"+" <- peakSignals[["+"]][idxbfPeakPass,]
      bfPassPeakSignals$"-" <- peakSignals[["-"]][idxbfPeakPass,]
      bfPassPeakSites <- peakSites[idxbfPeakPass]
      numbfPassPeakSites <- length(bfPassPeakSignals[["+"]][,1])
      bfPassPeakSignalTotals <- genomeSignalTotals[idxbfPeakPass]
      bfPassPeakMotifSignalTotals <- motifGenomeSignalTotals[idxbfPeakPass]
      
      ## Generate combined signal for bf passing sites
      combinedbfPassPeakSignal <- list()
      combinedbfPassPeakSignal$signal <- matrix(data = NA, nrow = numbfPassPeakSites, ncol = bpTotal)
      for (i in 1:numbfPassPeakSites){combinedbfPassPeakSignal$signal[i,] <- bfPassPeakSignals[["+"]][i,] + bfPassPeakSignals[["-"]][i,]}
      
      ## Generate profile for bf corrected peak sites
      cat("Generating profile for bf corrected peak sites", "\n")
      bfPeakProfile <- buildProfile(bfPassPeakSignals, libraryFactor, upstream, motifWidth, downstream)
      bfprof <- bfPeakProfile$profile
      
      ## Make graph for bf passing sites
      svg_path <- paste0(dirpath, "footprints/graphs/", samplename, ".", genename, ".", "motif", x, ".bfpeak.sites.svg")
      svg(file = svg_path) # set the filepath for saving the svg figure
      cat("Saving svg footprint image at path:", svg_path, "\n")
      PWMin <- pwm2pfm(PWM)
      cat("Plotting graph for bf corrected peaks...", "\n")
      ATACseqQC:::plotFootprints(bfprof,
                                 Mlen=motifWidth, motif=PWMin)
      dev.off()
      
      
      #### Make heatmap for bf passing sites ####
      ## USE THIS STRUCTURE FOR HEATMAPS
      ## first, combine the signals from plus and minus strand
      cat("Generating heatmap...", "\n")
      
      cat("2", "\n")
      heatsigsPlus <- bfPassPeakSignals$"+"
      heatsigsMinus <- bfPassPeakSignals$"-"

      cat("3", "\n")
      heatnumsites <- length(heatsigsPlus[,1])
      heatnumbp <- length(heatsigsPlus[1,])

      cat("4", "\n")
      cat(heatnumbp, "\n")
      cat(heatnumsites, "\n")
      
      cat("5", "\n")
      heatSites <- bfPassPeakSites
      
      ##
      temp <- matrix(data = NA, nrow = heatnumsites, ncol = heatnumbp)
      cat("6", "\n")
      
      for (i in 1:heatnumsites){temp[i,] <- heatsigsPlus[i,] + heatsigsMinus[i,]}
      cat("7", "\n")
      
      ## scale each row individually
      for (a in 1:heatnumsites){
        maxsig <- max(temp[a,])
        for (b in 1:heatnumbp){temp[a,b] <- (temp[a,b] / maxsig)}}
      cat("9", "\n")
      
      ## Set the maximal value of each datapoint to 1
      maxsig <- 1
      
      ## invert signals
      for (a in 1:heatnumsites){for (b in 1:heatnumbp){temp[a,b] <- (1-temp[a,b])}}
      cat("10", "\n")
      
      ## Annotate the combined sublist name which will become the tital of the heatmap plot
      cat("11", "\n")
      plottitle <- paste0(genename, "_motif", x, "_numsites", heatnumsites)
      
      combined <- list()
      com <- paste0("combined$", plottitle, " <- temp")
      eval(parse(text = com))
      
      ##
      cat("12", "\n")
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
      cat("13", "\n")
      ChIPpeakAnno::featureAlignedHeatmap(combined,
                                          feature.gr = reCenterPeaks(heatSites, width = heatnumbp),
                                          upper.extreme = maxsig, # set this to control the heatmap scale
                                          annoMcols = "score",
                                          sortBy = "score",
                                          n.tile = heatnumbp,
                                          margin = c(0.1, 0.005, 0.05, 0.2),
                                          color = colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias = 0.9)(100),
                                          gp = gpar(fontsize=10),
                                          newpage = TRUE)
      dev.off()
      ############################################################
      
      
      
    }, # end try
    error=function(cond){
      message(cond)
      return(NA)
    },
    finally={})
    
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
    parsedSitesInfo$heatsites <- heatSites
    parsedSitesInfo$heatnumbp <- heatnumbp
    ##
    save(parsedSitesInfo, file = info_path)

  } # end if (file.exists(info_path) == TRUE)
} # end for (x in 1:num_motifs)

cat("Finished...", "\n")
file.create(outpathdone)