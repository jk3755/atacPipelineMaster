
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)


#### Library loads
cat("Loading libraries...", "\n")
library(ATACseqQC)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)


#### Functions
cat("Creating functions...", "\n")
getMotifPWM <- function(symbol)
{
  
  #### Query MotifDb ####
  # Define string to return only Hsapiens motifs
  organism_rows = grep('Hsapiens', values(MotifDb)$organism, ignore.case = TRUE)
  # Define string for given gene
  gene_symbol_rows = grep(symbol, values(MotifDb)$geneSymbol, ignore.case = TRUE)
  # Get indices for the intersection of gene and organism
  human_gene_rows = intersect(gene_symbol_rows, organism_rows)
  # Pull the PWMs
  # need to make it a list to use unique() function to remove duplicate entries
  gene_motifs <- as.list(MotifDb[human_gene_rows])
  # Get unique motifs
  unique_gene_motifs <- unique(gene_motifs)
  # Get number of motifs
  number_unique_gene_motifs <- length(unique_gene_motifs)
  #### end ####
  
  #### Return PWM ####
  # Make a vector containing all the motifs
  unique_gene_motifs.vector <- c()
  for (a in 1:number_unique_gene_motifs) {
    unique_gene_motifs.vector[a] <- unique_gene_motifs[a]
  }
  
  # return the vector
  return(unique_gene_motifs.vector)
  #### end ####
  
} # end getMotifPWM function


#### Get the PWM for current gene
# pull the current gene symbol from snakemake
gene <- snakemake@wildcards[["gene"]]
cat("Current gene is: ", gene, "\n")
cat("Fetching PWMs...", "\n")
# get all human annotated records
PWMlist <- getMotifPWM(gene)


#### Set parameters
cat("Setting input parameters...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
motif_score <- "90%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens
PWM <- PWMlist[[1]]


#### Generate the footprint
cat("Generating footprint signal...", "\n")
# generate signal
signal <- factorFootprints(
   # bam file input
   bamfiles = bampath,
   # bai index input
   index = baipath,
   # PWM input
   pfm = PWM,
   # reference genome
   genome = genome,
   # minimum motif matching score
   min.score = motif_score,
   # scope of analysis
   seqlev = scope,
   # bp upstream of motif to analyze
   upstream = upstream,
   # bp downstream of motif to analyze
   downstream = downstream
)


#### Save
cat("Saving...", "\n")
outpath <- snakemake@output[[1]]
save(outpath)


## pseudocode ##
# 1) load sites object and use it to calculate wg signal
# 2) Find total signal at every site
# 3) Take top 10% sites
# 4) Find all unique values for total signal
# 5) Generate a null model for each unique total signal
# 6) Perform t-test and generate p-values
# 7) Get bonferroni corrected p-value passing sites
# 8) Subset binding sites and regenerate profile
# 9) Plot footprint and save profiles

#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)

#### Library loads
gc()
cat("Loading libraries...", "\n")
suppressMessages(library(MotifDb))
suppressMessages(library(ATACseqQC))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ChIPpeakAnno))

#### Snakemake variables
cat("Loading snakemake variables...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitesfile <- snakemake@input[[3]]
out_filepath <- snakemake@output[[1]]
sample_name <- snakemake@wildcards[["sample"]]
gene_name <- snakemake@wildcards[["gene"]]
cat("Sample name:", sample_name, "\n", "Gene name:", gene_name, "\n")

#### Load workspace
cat("Loading workspace...", "\n")
load(sitesfile)
num_motif <- length(bindingSites)
rm(bindingSites, genome, PWM, PWMlist, sites, temp, b, gene, num_pwm, outpath, getMotifPWM)
gc()

#### Parameters
cat("Setting input parameters...", "\n")
motif_score <- "95%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens

#### Functions
cat("Building functions...", "\n")
generateNullFP <- function(n, t, s, m){
  #This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
  #- Consider the total signal (number of insertions) at each specific ~200 bp locus
  #- Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
  #- Use published or experimentally derived models of Tn5 sequence specific insertion bias
  #- For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
  #- Generate the null model graph by weighted random residstribution of the total observed signal at that site
  #- Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
  #These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
  
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
#
makeProfile <- function(sites, PWM, type){
  
  anchor <- c("cut site")
  pfm <- PWM  
  min.score <- motif_score
  bindingSites <- sites
  seqlev <- scope
  bamfiles <- bampath
  index <- baipath
  
  stopifnot(is(bindingSites, "GRanges"))
  stopifnot(all(!is.na(seqlengths(bindingSites))))
  stopifnot(length(bindingSites)>1)
  stopifnot(length(bindingSites$score)==length(bindingSites))
  mt <- bindingSites
  mt$userdefined <- TRUE
  
  wid <- ncol(pfm)
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which)
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
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
  
  
  ## split into positive strand and negative strand
  bamIn <- split(bamIn, strand(bamIn))
  
  
  ## get coverage
  cvglist <- sapply(bamIn, coverage)
  cvglist <- cvglist[c("+", "-")]
  cvglist <- lapply(cvglist, function(.ele)
    .ele[names(.ele) %in% seqlev])
  
  ## coverage of mt, must be filtered, otherwise too much
  cvgSum <- cvglist[["+"]] + cvglist[["-"]]
  mt.s <- split(mt, seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  
  ## too much if use upstream and downstream, just use 3*wid maybe better.
  mt.s.ext <- promoters(mt.s, upstream=wid, downstream=wid+wid)
  stopifnot(all(lengths(mt.s.ext)==lengths(mt.s)))
  mt.v <- Views(cvgSum, mt.s.ext)
  mt.s <- mt.s[viewSums(mt.v)>0] 
  mt <- unlist(mt.s)
  mt.ids <- promoters(reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2)+1)
  mt.ids <- paste0(as.character(seqnames(mt.ids)), ":", start(mt.ids), "-", end(mt.ids))
  sigs <- featureAlignedSignal(cvglists=cvglist,
                               feature.gr=reCenterPeaks(mt, width=1),
                               upstream=upstream+floor(wid/2),
                               downstream=downstream+ceiling(wid/2),
                               n.tile=upstream+downstream+wid)
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  cor <- lapply(sigs, function(sig){
    sig.colMeans <- colMeans(sig)
    ## calculate correlation of footprinting and binding score
    windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, step = 1)[[1]]
    # remove the windows with overlaps of motif binding region
    windows <- windows[end(windows)<=upstream | start(windows)>=upstream+wid]
    sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
    windows.sel <- windows[which.max(sig.windowMeans)][1]
    highest.sig.windows <- 
      rowMeans(sig[, start(windows.sel):end(windows.sel)])
    predictedBindingSiteScore <- mt$score
    if(length(predictedBindingSiteScore) == length(highest.sig.windows)){
      suppressWarnings({
        cor <- cor.test(x = predictedBindingSiteScore, 
                        y = highest.sig.windows, 
                        method = "spearman")
      })
    }else{
      cor <- NA
    }
    cor
  })
  
  ###########################################################
  ## This is where you can insert the shuffled signals ######
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  ## Make the null model
  if (type == "nullmodel"){
    nums <- nrow(sigs[["+"]])
    cols <- ncol(sigs[["+"]])
    plusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    minusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    for (c in 1:nums){
      plusshuf[c,] <- as.double(sample(sigs[["+"]][c,], replace=FALSE))
      minusshuf[c,] <- as.double(sample(sigs[["-"]][c,], replace=FALSE))
    }
    sigs$"+" <- plusshuf
    sigs$"-" <- minusshuf
  }
  ###########################################################
  ###########################################################
  
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  
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
  
  pwm2pfm <- function(pfm, name="motif"){
    if(!all(round(colSums(pfm), digits=4)==1)){
      return(NULL)
    }
    new("pfm", mat=as.matrix(pfm), name=name)
  }
  
  ## Input arguments for calling the plotFootprints function
  args <- list()
  args$Profile <- c(Profile[["+"]], Profile[["-"]])
  args$ylab <- ifelse(anchor=="cut site", "Cut-site probability", "reads density (arbitrary unit)")
  args$Mlen <- wid
  args$motif <- pwm2pfm(pfm)
  args$segmentation <- Profile.seg
  
  if (type == "nullmodel"){
    args$nullplus <- sigs$"+"
    args$nullminus <- sigs$"-"
  }
  
  return(args)
}
#
pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}

#### Generate WG signals with binding sites object
for (a in 1:num_motif){ # An iterator for each unique motif for this gene
  
  #
  cat("Retrieving stored binding sites...", "\n")
  load(sitesfile)
  PWM <- bindingSites[[a]][["PWM"]]
  sites <- bindingSites[[a]][["sites"]]
  rm(bindingSites, PWMlist, temp, b, gene, num_pwm, outpath, getMotifPWM)
  gc()
  
  #
  cat("Generating WG footprint signal", a, "of ", num_motif,  "\n")
  sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
  sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
  
  # generate signal
  wgsignals <- factorFootprints(
    # bam file input
    bamfiles = bampath,
    # bai index input
    index = baipath,
    # PWM input
    pfm = PWM,
    # reference genome
    genome = genome,
    # minimum motif matching score
    min.score = motif_score,
    # scope of analysis
    seqlev = scope,
    # bp upstream of motif to analyze
    upstream = upstream,
    # bp downstream of motif to analyze
    downstream = downstream,
    # binding sites
    bindingSites = sites
  )
  
  # get the number of sites in the whole set
  total_sites <- nrow(wgsignals[["signal"]][["+"]])
  cat("Total sites found:", total_sites, "\n")
  
  ## sum up the signals
  ## NOTE - preserve the row names = genomic coordinates
  plus <- wgsignals[["signal"]][["+"]]
  minus <-  wgsignals[["signal"]][["-"]]
  # combine the plus and minus strand signals
  total_bp <- length(plus[1,])
  combined <- matrix(data=NA,nrow=total_sites,ncol=total_bp)
  # combine the signals
  for (i in 1:total_sites){combined[i,] <- plus[i,] + minus[i,]}
  # xfer row names
  rownames(combined) <- rownames(plus)
  # cleanup plus and minus
  rm(plus,minus)
  ## get the motif width
  motif_width <- total_bp - 200
  
  ## for each site, calculate the total signal and motif signal
  cat("Calculating signal at each site...", "\n")
  signal_totals <- c() # total signal
  motifsignal_totals <- c() # total signal in the motif
  for (i in 1:total_sites){
    signal_totals[i] <- sum(combined[i,])
    motifsignal_totals[i] <- sum(combined[i,100:(100+motif_width)])}
  
  ## perform quantile analysis and get top 10%
  # retrieve the specified quantiles from the total signals
  quantiles <- quantile(signal_totals, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  # retrieve the indices of the sites in the top 10% of total signal
  topten_ind <- which(signal_totals > quantiles[10])
  # subset the sites data
  topten_sites <- combined[topten_ind,]
  # xfer the rownames
  rownames(topten_sites) <- rownames(combined[topten_ind,])
  # reset the number of inputs
  topten_numsites <- nrow(topten_sites)
  ## recalculate total signal at each site
  topten_signal_totals <- signal_totals[topten_ind] # total signal
  topten_motifsignal_totals <- motifsignal_totals[topten_ind] # total signal in the motif
  
  ## get all the unique values for total signal in the threshold passing sites
  ## this will be used to compute the null models
  unique_signals <- unique(topten_signal_totals)
  
  ## generate the set of null models
  ## Remember, you must store the null models for proper record keeping and reproducibility
  cat("Generating null models...", "\n")
  null_sim <- list()
  num <- length(unique_signals)
  null_mean <- c()
  for (i in 1:num){
    n <- 1000
    t <- unique_signals[i]
    s <- total_bp
    m <- motif_width
    null_sim[[i]] <- generateNullFP(n,t,s,m)
    null_mean[i] <- mean(null_sim[[i]])}
  
  ## Perform a one-tailed t-test to generate a p-value for each observed motif site
  ttest <- list() # list to store the results of the t-tests
  pvalue <- c() # vector to store the p-values
  tvalue <- c() # vector to store the t-value
  
  for (i in 1:topten_numsites){
    current <- c(topten_sites[i,])
    signal <- c(sum(current))
    # retrieve the appropriate null model
    nullmodel <- null_sim[[which(unique_signals==signal)]]
    # do the t-test
    ttest[[i]] <- t.test(current[100:116], mu=mean(nullmodel), alternative="less", conf.level = 0.95)
    pvalue[i] <- ttest[[i]][["p.value"]]
    tvalue[i] <- ttest[[i]][["statistic"]][["t"]]}
  
  # get the indices of the sites that are lower than p = 0.05
  ppass_ind <- which(pvalue < 0.05)
  pvalue_pass <- pvalue[ppass_ind]
  
  # bonferroni correction
  bf_ppass_ind <- which(pvalue < (0.05/topten_numsites))
  bf_pvalue_pass <- pvalue[bf_ppass_ind]
  
  ## Subset the wg signals object bindingSites
  ## make profile function only requires binding sites and PWM as input
  #pvalue_pass_bindingsites <- topten_bindingsites[ppass_ind]
  #bh_pass_bindingsites <- topten_bindingsites[bh_pass_ind]
  topten_bindingsites <- wgsignals[["bindingSites"]][topten_ind]
  bf_pass_bindingsites <- topten_bindingsites[bf_ppass_ind]
  
  # generate the profiles for all objects
  cat("Generating CENTIPEDE profile...", "\n")
  #topten_profile <- makeProfile(topten_bindingsites, PWM, type="normal")
  #pvalue_profile <- makeProfile(pvalue_pass_bindingsites, PWM, type="normal")
  bf_profile <- makeProfile(bf_pass_bindingsites, PWM, type="normal")
  #bh_profile <- makeProfile(bh_pass_bindingsites, PWM, type="normal")
  
  #
  svg_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".motif", a, ".svg")
  svg(file = svg_path) # set the filepath for saving the svg figure
  cat("Saving svg footprint image at path:", svg_path, "\n")
  
  ## FP graph
  PWMin <- pwm2pfm(PWM)
  cat("Plotting graph...", "\n")
  ATACseqQC:::plotFootprints(bf_profile$Profile,
                             Mlen=motif_width, motif=PWMin)
  dev.off()
  
  ##
  cat("Transferring data to FPdata object...", "\n")
  FPdata <- list()
  FPdata$motif <- a
  FPdata$PWM <- PWM
  FPdata$wg_bindingsites <- sites
  FPdata$wg_numsites <- total_sites
  FPdata$wg_signals <- wgsignals
  FPdata$bf_pass_bindingsites <- bf_pass_bindingsites
  FPdata$bf_pvalue_pass <- bf_pvalue_pass
  FPdata$bf_profile <- bf_profile
  FPdata$pvalue <- pvalue
  FPdata$tvalue <- tvalue
  FPdata$null_sim <- null_sim
  FPdata$wg_top10_signaltotals <- topten_signal_totals
  FPdata$wg_top10_motifsignal_totals <- topten_motifsignal_totals
  FPdata$wg_signaltotals <- signal_totals
  FPdata$wg_motifsignal_totals <- motifsignal_totals
  ##
  rm(PWM, sites, total_sites, wgsignals, bf_pass_bindingsites, bf_pvalue_pass, bf_profile, pvalue, tvalue, null_sim, topten_signal_totals, topten_motifsignal_totals, signal_totals, motifsignal_totals)
  signal_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".", "motif", a, ".Rdata")
  cat("Saving signals data for gene", gene_name, "motif", a, "at path:", signal_path, "\n")
  save.image(file = signal_path)
  rm(FPdata)
  gc()
  
}

############

#### Save
cat("Finalizing...\n")
file.create(file = out_filepath)

#### Cleanup workspace
cat("Cleaning workspace...", "\n")
rm(list=ls())
gc()

## pseudocode ##
# 1) load sites object and use it to calculate wg signal
# 2) Find total signal at every site
# 3) Take top 10% sites
# 4) Find all unique values for total signal
# 5) Generate a null model for each unique total signal
# 6) Perform t-test and generate p-values
# 7) Get bonferroni corrected p-value passing sites
# 8) Subset binding sites and regenerate profile
# 9) Plot footprint and save profiles

#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)

#### Library loads
gc()
cat("Loading libraries...", "\n")
suppressMessages(library(MotifDb))
suppressMessages(library(ATACseqQC))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ChIPpeakAnno))

#### Snakemake variables
cat("Loading snakemake variables...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitesfile <- snakemake@input[[3]]
out_filepath <- snakemake@output[[1]]
sample_name <- snakemake@wildcards[["sample"]]
gene_name <- snakemake@wildcards[["gene"]]
cat("Sample name:", sample_name, "\n", "Gene name:", gene_name, "\n")

#### Load workspace
cat("Loading workspace...", "\n")
load(sitesfile)
num_motif <- length(bindingSites)
rm(bindingSites, genome, PWM, PWMlist, sites, temp, b, gene, num_pwm, outpath, getMotifPWM)
gc()

#### Parameters
cat("Setting input parameters...", "\n")
motif_score <- "95%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens

#### Functions
cat("Building functions...", "\n")
generateNullFP <- function(n, t, s, m){
  #This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
  #- Consider the total signal (number of insertions) at each specific ~200 bp locus
  #- Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
  #- Use published or experimentally derived models of Tn5 sequence specific insertion bias
  #- For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
  #- Generate the null model graph by weighted random residstribution of the total observed signal at that site
  #- Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
  #These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
  
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
#
makeProfile <- function(sites, PWM, type){
  
  anchor <- c("cut site")
  pfm <- PWM  
  min.score <- motif_score
  bindingSites <- sites
  seqlev <- scope
  bamfiles <- bampath
  index <- baipath
  
  stopifnot(is(bindingSites, "GRanges"))
  stopifnot(all(!is.na(seqlengths(bindingSites))))
  stopifnot(length(bindingSites)>1)
  stopifnot(length(bindingSites$score)==length(bindingSites))
  mt <- bindingSites
  mt$userdefined <- TRUE
  
  wid <- ncol(pfm)
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which)
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
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
  
  
  ## split into positive strand and negative strand
  bamIn <- split(bamIn, strand(bamIn))
  
  
  ## get coverage
  cvglist <- sapply(bamIn, coverage)
  cvglist <- cvglist[c("+", "-")]
  cvglist <- lapply(cvglist, function(.ele)
    .ele[names(.ele) %in% seqlev])
  
  ## coverage of mt, must be filtered, otherwise too much
  cvgSum <- cvglist[["+"]] + cvglist[["-"]]
  mt.s <- split(mt, seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  
  ## too much if use upstream and downstream, just use 3*wid maybe better.
  mt.s.ext <- promoters(mt.s, upstream=wid, downstream=wid+wid)
  stopifnot(all(lengths(mt.s.ext)==lengths(mt.s)))
  mt.v <- Views(cvgSum, mt.s.ext)
  mt.s <- mt.s[viewSums(mt.v)>0] 
  mt <- unlist(mt.s)
  mt.ids <- promoters(reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2)+1)
  mt.ids <- paste0(as.character(seqnames(mt.ids)), ":", start(mt.ids), "-", end(mt.ids))
  sigs <- featureAlignedSignal(cvglists=cvglist,
                               feature.gr=reCenterPeaks(mt, width=1),
                               upstream=upstream+floor(wid/2),
                               downstream=downstream+ceiling(wid/2),
                               n.tile=upstream+downstream+wid)
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  cor <- lapply(sigs, function(sig){
    sig.colMeans <- colMeans(sig)
    ## calculate correlation of footprinting and binding score
    windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, step = 1)[[1]]
    # remove the windows with overlaps of motif binding region
    windows <- windows[end(windows)<=upstream | start(windows)>=upstream+wid]
    sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
    windows.sel <- windows[which.max(sig.windowMeans)][1]
    highest.sig.windows <- 
      rowMeans(sig[, start(windows.sel):end(windows.sel)])
    predictedBindingSiteScore <- mt$score
    if(length(predictedBindingSiteScore) == length(highest.sig.windows)){
      suppressWarnings({
        cor <- cor.test(x = predictedBindingSiteScore, 
                        y = highest.sig.windows, 
                        method = "spearman")
      })
    }else{
      cor <- NA
    }
    cor
  })
  
  ###########################################################
  ## This is where you can insert the shuffled signals ######
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  ## Make the null model
  if (type == "nullmodel"){
    nums <- nrow(sigs[["+"]])
    cols <- ncol(sigs[["+"]])
    plusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    minusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    for (c in 1:nums){
      plusshuf[c,] <- as.double(sample(sigs[["+"]][c,], replace=FALSE))
      minusshuf[c,] <- as.double(sample(sigs[["-"]][c,], replace=FALSE))
    }
    sigs$"+" <- plusshuf
    sigs$"-" <- minusshuf
  }
  ###########################################################
  ###########################################################
  
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  
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
  
  pwm2pfm <- function(pfm, name="motif"){
    if(!all(round(colSums(pfm), digits=4)==1)){
      return(NULL)
    }
    new("pfm", mat=as.matrix(pfm), name=name)
  }
  
  ## Input arguments for calling the plotFootprints function
  args <- list()
  args$Profile <- c(Profile[["+"]], Profile[["-"]])
  args$ylab <- ifelse(anchor=="cut site", "Cut-site probability", "reads density (arbitrary unit)")
  args$Mlen <- wid
  args$motif <- pwm2pfm(pfm)
  args$segmentation <- Profile.seg
  
  if (type == "nullmodel"){
    args$nullplus <- sigs$"+"
    args$nullminus <- sigs$"-"
  }
  
  return(args)
}
#
pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}

#### Generate WG signals with binding sites object
for (a in 1:num_motif){ # An iterator for each unique motif for this gene
  
  signal_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".", "motif", a, ".Rdata")
  if (file.exists(signal_path) == TRUE){
    cat("File already exists, skipping...", "\n")
    next
  } else {
    
    #
    cat("Retrieving stored binding sites...", "\n")
    load(sitesfile)
    PWM <- bindingSites[[a]][["PWM"]]
    sites <- bindingSites[[a]][["sites"]]
    rm(bindingSites, PWMlist, temp, b, gene, num_pwm, outpath, getMotifPWM)
    gc()
    
    #
    cat("Generating WG footprint signal", a, "of ", num_motif,  "\n")
    sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
    sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
    
    # generate signal
    wgsignals <- factorFootprints(
      # bam file input
      bamfiles = bampath,
      # bai index input
      index = baipath,
      # PWM input
      pfm = PWM,
      # reference genome
      genome = genome,
      # minimum motif matching score
      min.score = motif_score,
      # scope of analysis
      seqlev = scope,
      # bp upstream of motif to analyze
      upstream = upstream,
      # bp downstream of motif to analyze
      downstream = downstream,
      # binding sites
      bindingSites = sites
    )
    
    # get the number of sites in the whole set
    total_sites <- nrow(wgsignals[["signal"]][["+"]])
    cat("Total sites found:", total_sites, "\n")
    
    ## sum up the signals
    ## NOTE - preserve the row names = genomic coordinates
    plus <- wgsignals[["signal"]][["+"]]
    minus <-  wgsignals[["signal"]][["-"]]
    # combine the plus and minus strand signals
    total_bp <- length(plus[1,])
    combined <- matrix(data=NA,nrow=total_sites,ncol=total_bp)
    # combine the signals
    for (i in 1:total_sites){combined[i,] <- plus[i,] + minus[i,]}
    # xfer row names
    rownames(combined) <- rownames(plus)
    # cleanup plus and minus
    rm(plus,minus)
    ## get the motif width
    motif_width <- total_bp - 200
    
    ## for each site, calculate the total signal and motif signal
    cat("Calculating signal at each site...", "\n")
    signal_totals <- c() # total signal
    motifsignal_totals <- c() # total signal in the motif
    for (i in 1:total_sites){
      signal_totals[i] <- sum(combined[i,])
      motifsignal_totals[i] <- sum(combined[i,100:(100+motif_width)])}
    
    ## perform quantile analysis and get top 10%
    # retrieve the specified quantiles from the total signals
    quantiles <- quantile(signal_totals, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
    # retrieve the indices of the sites in the top 10% of total signal
    topten_ind <- which(signal_totals > quantiles[10])
    # subset the sites data
    topten_sites <- combined[topten_ind,]
    # xfer the rownames
    rownames(topten_sites) <- rownames(combined[topten_ind,])
    # reset the number of inputs
    topten_numsites <- nrow(topten_sites)
    ## recalculate total signal at each site
    topten_signal_totals <- signal_totals[topten_ind] # total signal
    topten_motifsignal_totals <- motifsignal_totals[topten_ind] # total signal in the motif
    
    ## get all the unique values for total signal in the threshold passing sites
    ## this will be used to compute the null models
    unique_signals <- unique(topten_signal_totals)
    
    ## generate the set of null models
    ## Remember, you must store the null models for proper record keeping and reproducibility
    cat("Generating null models...", "\n")
    null_sim <- list()
    num <- length(unique_signals)
    null_mean <- c()
    for (i in 1:num){
      n <- 1000
      t <- unique_signals[i]
      s <- total_bp
      m <- motif_width
      null_sim[[i]] <- generateNullFP(n,t,s,m)
      null_mean[i] <- mean(null_sim[[i]])}
    
    ## Perform a one-tailed t-test to generate a p-value for each observed motif site
    ttest <- list() # list to store the results of the t-tests
    pvalue <- c() # vector to store the p-values
    tvalue <- c() # vector to store the t-value
    
    for (i in 1:topten_numsites){
      current <- c(topten_sites[i,])
      signal <- c(sum(current))
      # retrieve the appropriate null model
      nullmodel <- null_sim[[which(unique_signals==signal)]]
      # do the t-test
      ttest[[i]] <- t.test(current[100:116], mu=mean(nullmodel), alternative="less", conf.level = 0.95)
      pvalue[i] <- ttest[[i]][["p.value"]]
      tvalue[i] <- ttest[[i]][["statistic"]][["t"]]}
    
    # get the indices of the sites that are lower than p = 0.05
    ppass_ind <- which(pvalue < 0.05)
    pvalue_pass <- pvalue[ppass_ind]
    
    # bonferroni correction
    bf_ppass_ind <- which(pvalue < (0.05/topten_numsites))
    bf_pvalue_pass <- pvalue[bf_ppass_ind]
    
    ## Subset the wg signals object bindingSites
    ## make profile function only requires binding sites and PWM as input
    #pvalue_pass_bindingsites <- topten_bindingsites[ppass_ind]
    #bh_pass_bindingsites <- topten_bindingsites[bh_pass_ind]
    topten_bindingsites <- wgsignals[["bindingSites"]][topten_ind]
    bf_pass_bindingsites <- topten_bindingsites[bf_ppass_ind]
    
    # generate the profiles for all objects
    cat("Generating CENTIPEDE profile...", "\n")
    #topten_profile <- makeProfile(topten_bindingsites, PWM, type="normal")
    #pvalue_profile <- makeProfile(pvalue_pass_bindingsites, PWM, type="normal")
    bf_profile <- makeProfile(bf_pass_bindingsites, PWM, type="normal")
    #bh_profile <- makeProfile(bh_pass_bindingsites, PWM, type="normal")
    
    #
    svg_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".motif", a, ".svg")
    svg(file = svg_path) # set the filepath for saving the svg figure
    cat("Saving svg footprint image at path:", svg_path, "\n")
    
    ## FP graph
    PWMin <- pwm2pfm(PWM)
    cat("Plotting graph...", "\n")
    ATACseqQC:::plotFootprints(bf_profile$Profile,
                               Mlen=motif_width, motif=PWMin)
    dev.off()
    
    ##
    cat("Transferring data to FPdata object...", "\n")
    FPdata <- list()
    FPdata$motif <- a
    FPdata$PWM <- PWM
    FPdata$wg_bindingsites <- sites
    FPdata$wg_numsites <- total_sites
    FPdata$wg_signals <- wgsignals
    FPdata$bf_pass_bindingsites <- bf_pass_bindingsites
    FPdata$bf_pvalue_pass <- bf_pvalue_pass
    FPdata$bf_profile <- bf_profile
    FPdata$pvalue <- pvalue
    FPdata$tvalue <- tvalue
    FPdata$null_sim <- null_sim
    FPdata$wg_top10_signaltotals <- topten_signal_totals
    FPdata$wg_top10_motifsignal_totals <- topten_motifsignal_totals
    FPdata$wg_signaltotals <- signal_totals
    FPdata$wg_motifsignal_totals <- motifsignal_totals
    ##
    rm(PWM, sites, total_sites, wgsignals, bf_pass_bindingsites, bf_pvalue_pass, bf_profile, pvalue, tvalue, null_sim, topten_signal_totals, topten_motifsignal_totals, signal_totals, motifsignal_totals)
    cat("Saving signals data for gene", gene_name, "motif", a, "at path:", signal_path, "\n")
    save.image(file = signal_path)
    rm(FPdata)
    gc()
  }
}

############

#### Save
cat("Finalizing...\n")
file.create(file = out_filepath)

#### Cleanup workspace
cat("Cleaning workspace...", "\n")
rm(list=ls())
gc()


## pseudocode ##
# 1) load sites object and use it to calculate wg signal
# 2) Find total signal at every site
# 3) Take top 10% sites
# 4) Find all unique values for total signal
# 5) Generate a null model for each unique total signal
# 6) Perform t-test and generate p-values
# 7) Get bonferroni corrected p-value passing sites
# 8) Subset binding sites and regenerate profile
# 9) Plot footprint and save profiles

#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)
#biocLite("Rsamtools", suppressUpdates = TRUE)

#### Library loads
gc()
cat("Loading libraries...", "\n")
suppressMessages(library(MotifDb))
suppressMessages(library(ATACseqQC))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(ChIPpeakAnno))

#### Snakemake variables
cat("Loading snakemake variables...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitesfile <- snakemake@input[[3]]
out_filepath <- snakemake@output[[1]]
sample_name <- snakemake@wildcards[["sample"]]
gene_name <- snakemake@wildcards[["gene"]]
cat("Sample name:", sample_name, "\n", "Gene name:", gene_name, "\n")

#### Load workspace
cat("Loading workspace...", "\n")
load(sitesfile)
num_motif <- length(bindingSites)
rm(bindingSites, genome, PWM, PWMlist, sites, temp, b, gene, num_pwm, outpath, getMotifPWM)
gc()

#### Parameters
cat("Setting input parameters...", "\n")
motif_score <- "95%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens

#### Functions
cat("Building functions...", "\n")
generateNullFP <- function(n, t, s, m){
  #This script will be used to generate indiviudal null models at predicted motif binding sites across the genome when scanning for TF footprinting from ATAC-seq data. To generate these null models, the current model will need to:
  #- Consider the total signal (number of insertions) at each specific ~200 bp locus
  #- Use the actul underlying reference sequence of that ~200 bp stretch from the hg38 reference genome
  #- Use published or experimentally derived models of Tn5 sequence specific insertion bias
  #- For each locus, build a probablistic model of insertion site distributions based on the underlying sequence and Tn5 insertion bias
  #- Generate the null model graph by weighted random residstribution of the total observed signal at that site
  #- Importantly, the null model must be generated separately for the plus and minus strand, it can then be combined and compared to the combined signal from the reference observed signal at that sequence
  #These null models can then be used for a site-by-site comparison of the null model against the observed data to accept or reject the null hypothesis
  
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
#
makeProfile <- function(sites, PWM, type){
  
  anchor <- c("cut site")
  pfm <- PWM  
  min.score <- motif_score
  bindingSites <- sites
  seqlev <- scope
  bamfiles <- bampath
  index <- baipath
  
  stopifnot(is(bindingSites, "GRanges"))
  stopifnot(all(!is.na(seqlengths(bindingSites))))
  stopifnot(length(bindingSites)>1)
  stopifnot(length(bindingSites$score)==length(bindingSites))
  mt <- bindingSites
  mt$userdefined <- TRUE
  
  wid <- ncol(pfm)
  seqlevels(mt) <- seqlev
  seqinfo(mt) <- Seqinfo(seqlev, seqlengths = seqlengths(mt))
  
  ## read in bam file with input seqlev specified by users
  which <- as(seqinfo(mt), "GRanges")
  param <- ScanBamParam(which=which)
  if(anchor=="cut site"){
    bamIn <- mapply(function(.b, .i) readGAlignments(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
  }else{
    bamIn <- mapply(function(.b, .i) readGAlignmentPairs(.b, .i, param = param), 
                    bamfiles, index, SIMPLIFY = FALSE)
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
  
  
  ## split into positive strand and negative strand
  bamIn <- split(bamIn, strand(bamIn))
  
  
  ## get coverage
  cvglist <- sapply(bamIn, coverage)
  cvglist <- cvglist[c("+", "-")]
  cvglist <- lapply(cvglist, function(.ele)
    .ele[names(.ele) %in% seqlev])
  
  ## coverage of mt, must be filtered, otherwise too much
  cvgSum <- cvglist[["+"]] + cvglist[["-"]]
  mt.s <- split(mt, seqnames(mt))
  seqlev <- intersect(names(cvgSum), names(mt.s))
  cvgSum <- cvgSum[seqlev]
  mt.s <- mt.s[seqlev]
  
  ## too much if use upstream and downstream, just use 3*wid maybe better.
  mt.s.ext <- promoters(mt.s, upstream=wid, downstream=wid+wid)
  stopifnot(all(lengths(mt.s.ext)==lengths(mt.s)))
  mt.v <- Views(cvgSum, mt.s.ext)
  mt.s <- mt.s[viewSums(mt.v)>0] 
  mt <- unlist(mt.s)
  mt.ids <- promoters(reCenterPeaks(mt, width=1),
                      upstream=upstream+floor(wid/2),
                      downstream=downstream+ceiling(wid/2)+1)
  mt.ids <- paste0(as.character(seqnames(mt.ids)), ":", start(mt.ids), "-", end(mt.ids))
  sigs <- featureAlignedSignal(cvglists=cvglist,
                               feature.gr=reCenterPeaks(mt, width=1),
                               upstream=upstream+floor(wid/2),
                               downstream=downstream+ceiling(wid/2),
                               n.tile=upstream+downstream+wid)
  mt <- mt[match(rownames(sigs[[1]]), mt.ids)]
  cor <- lapply(sigs, function(sig){
    sig.colMeans <- colMeans(sig)
    ## calculate correlation of footprinting and binding score
    windows <- slidingWindows(IRanges(1, ncol(sig)), width = wid, step = 1)[[1]]
    # remove the windows with overlaps of motif binding region
    windows <- windows[end(windows)<=upstream | start(windows)>=upstream+wid]
    sig.windowMeans <- viewMeans(Views(sig.colMeans, windows))
    windows.sel <- windows[which.max(sig.windowMeans)][1]
    highest.sig.windows <- 
      rowMeans(sig[, start(windows.sel):end(windows.sel)])
    predictedBindingSiteScore <- mt$score
    if(length(predictedBindingSiteScore) == length(highest.sig.windows)){
      suppressWarnings({
        cor <- cor.test(x = predictedBindingSiteScore, 
                        y = highest.sig.windows, 
                        method = "spearman")
      })
    }else{
      cor <- NA
    }
    cor
  })
  
  ###########################################################
  ## This is where you can insert the shuffled signals ######
  sigs <- lapply(sigs, function(.ele) .ele[mt$userdefined, ])
  ## Make the null model
  if (type == "nullmodel"){
    nums <- nrow(sigs[["+"]])
    cols <- ncol(sigs[["+"]])
    plusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    minusshuf <- matrix(data=NA, nrow=nums, ncol=cols)
    for (c in 1:nums){
      plusshuf[c,] <- as.double(sample(sigs[["+"]][c,], replace=FALSE))
      minusshuf[c,] <- as.double(sample(sigs[["-"]][c,], replace=FALSE))
    }
    sigs$"+" <- plusshuf
    sigs$"-" <- minusshuf
  }
  ###########################################################
  ###########################################################
  
  mt <- mt[mt$userdefined]
  mt$userdefined <- NULL
  
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
  
  pwm2pfm <- function(pfm, name="motif"){
    if(!all(round(colSums(pfm), digits=4)==1)){
      return(NULL)
    }
    new("pfm", mat=as.matrix(pfm), name=name)
  }
  
  ## Input arguments for calling the plotFootprints function
  args <- list()
  args$Profile <- c(Profile[["+"]], Profile[["-"]])
  args$ylab <- ifelse(anchor=="cut site", "Cut-site probability", "reads density (arbitrary unit)")
  args$Mlen <- wid
  args$motif <- pwm2pfm(pfm)
  args$segmentation <- Profile.seg
  
  if (type == "nullmodel"){
    args$nullplus <- sigs$"+"
    args$nullminus <- sigs$"-"
  }
  
  return(args)
}
#
pwm2pfm <- function(pfm, name="motif"){
  if(!all(round(colSums(pfm), digits=4)==1)){
    return(NULL)
  }
  new("pfm", mat=as.matrix(pfm), name=name)
}

#### Generate WG signals with binding sites object
for (a in 1:num_motif){ # An iterator for each unique motif for this gene
  
  signal_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".", "motif", a, ".Rdata")
  if (file.exists(signal_path) == TRUE){
    cat("File already exists, skipping...", "\n")
    next
  } else {
    
    #
    cat("Retrieving stored binding sites...", "\n")
    load(sitesfile)
    PWM <- bindingSites[[a]][["PWM"]]
    sites <- bindingSites[[a]][["sites"]]
    rm(bindingSites, PWMlist, temp, b, gene, num_pwm, outpath, getMotifPWM)
    gc()
    
    #
    cat("Generating WG footprint signal", a, "of ", num_motif,  "\n")
    sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
    sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
    
    # generate signal
    wgsignals <- factorFootprints(
      # bam file input
      bamfiles = bampath,
      # bai index input
      index = baipath,
      # PWM input
      pfm = PWM,
      # reference genome
      genome = genome,
      # minimum motif matching score
      min.score = motif_score,
      # scope of analysis
      seqlev = scope,
      # bp upstream of motif to analyze
      upstream = upstream,
      # bp downstream of motif to analyze
      downstream = downstream,
      # binding sites
      bindingSites = sites
    )
    
    # get the number of sites in the whole set
    total_sites <- nrow(wgsignals[["signal"]][["+"]])
    cat("Total sites found:", total_sites, "\n")
    
    ## sum up the signals
    ## NOTE - preserve the row names = genomic coordinates
    plus <- wgsignals[["signal"]][["+"]]
    minus <-  wgsignals[["signal"]][["-"]]
    # combine the plus and minus strand signals
    total_bp <- length(plus[1,])
    combined <- matrix(data=NA,nrow=total_sites,ncol=total_bp)
    # combine the signals
    for (i in 1:total_sites){combined[i,] <- plus[i,] + minus[i,]}
    # xfer row names
    rownames(combined) <- rownames(plus)
    # cleanup plus and minus
    rm(plus,minus)
    ## get the motif width
    motif_width <- total_bp - 200
    
    ## for each site, calculate the total signal and motif signal
    cat("Calculating signal at each site...", "\n")
    signal_totals <- c() # total signal
    motifsignal_totals <- c() # total signal in the motif
    for (i in 1:total_sites){
      signal_totals[i] <- sum(combined[i,])
      motifsignal_totals[i] <- sum(combined[i,100:(100+motif_width)])}
    
    ## perform quantile analysis and get top 10%
    # retrieve the specified quantiles from the total signals
    quantiles <- quantile(signal_totals, probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
    # retrieve the indices of the sites in the top 10% of total signal
    topten_ind <- which(signal_totals > quantiles[10])
    # subset the sites data
    topten_sites <- combined[topten_ind,]
    # xfer the rownames
    rownames(topten_sites) <- rownames(combined[topten_ind,])
    # reset the number of inputs
    topten_numsites <- nrow(topten_sites)
    ## recalculate total signal at each site
    topten_signal_totals <- signal_totals[topten_ind] # total signal
    topten_motifsignal_totals <- motifsignal_totals[topten_ind] # total signal in the motif
    
    ## get all the unique values for total signal in the threshold passing sites
    ## this will be used to compute the null models
    unique_signals <- unique(topten_signal_totals)
    
    ## generate the set of null models
    ## Remember, you must store the null models for proper record keeping and reproducibility
    cat("Generating null models...", "\n")
    null_sim <- list()
    num <- length(unique_signals)
    null_mean <- c()
    for (i in 1:num){
      n <- 1000
      t <- unique_signals[i]
      s <- total_bp
      m <- motif_width
      null_sim[[i]] <- generateNullFP(n,t,s,m)
      null_mean[i] <- mean(null_sim[[i]])}
    
    ## Perform a one-tailed t-test to generate a p-value for each observed motif site
    ttest <- list() # list to store the results of the t-tests
    pvalue <- c() # vector to store the p-values
    tvalue <- c() # vector to store the t-value
    
    for (i in 1:topten_numsites){
      current <- c(topten_sites[i,])
      signal <- c(sum(current))
      # retrieve the appropriate null model
      nullmodel <- null_sim[[which(unique_signals==signal)]]
      # do the t-test
      ttest[[i]] <- t.test(current[100:116], mu=mean(nullmodel), alternative="less", conf.level = 0.95)
      pvalue[i] <- ttest[[i]][["p.value"]]
      tvalue[i] <- ttest[[i]][["statistic"]][["t"]]}
    
    # get the indices of the sites that are lower than p = 0.05
    ppass_ind <- which(pvalue < 0.05)
    pvalue_pass <- pvalue[ppass_ind]
    
    # bonferroni correction
    bf_ppass_ind <- which(pvalue < (0.05/topten_numsites))
    bf_pvalue_pass <- pvalue[bf_ppass_ind]
    
    ## Subset the wg signals object bindingSites
    ## make profile function only requires binding sites and PWM as input
    #pvalue_pass_bindingsites <- topten_bindingsites[ppass_ind]
    #bh_pass_bindingsites <- topten_bindingsites[bh_pass_ind]
    topten_bindingsites <- wgsignals[["bindingSites"]][topten_ind]
    bf_pass_bindingsites <- topten_bindingsites[bf_ppass_ind]
    
    # generate the profiles for all objects
    cat("Generating CENTIPEDE profile...", "\n")
    #topten_profile <- makeProfile(topten_bindingsites, PWM, type="normal")
    #pvalue_profile <- makeProfile(pvalue_pass_bindingsites, PWM, type="normal")
    bf_profile <- makeProfile(bf_pass_bindingsites, PWM, type="normal")
    #bh_profile <- makeProfile(bh_pass_bindingsites, PWM, type="normal")
    
    #
    svg_path <- paste0("coad_footprints/", sample_name, ".", gene_name, ".motif", a, ".svg")
    svg(file = svg_path) # set the filepath for saving the svg figure
    cat("Saving svg footprint image at path:", svg_path, "\n")
    
    ## FP graph
    PWMin <- pwm2pfm(PWM)
    cat("Plotting graph...", "\n")
    ATACseqQC:::plotFootprints(bf_profile$Profile,
                               Mlen=motif_width, motif=PWMin)
    dev.off()
    
    ##
    cat("Transferring data to FPdata object...", "\n")
    FPdata <- list()
    FPdata$motif <- a
    FPdata$PWM <- PWM
    FPdata$wg_bindingsites <- sites
    FPdata$wg_numsites <- total_sites
    FPdata$wg_signals <- wgsignals
    FPdata$bf_pass_bindingsites <- bf_pass_bindingsites
    FPdata$bf_pvalue_pass <- bf_pvalue_pass
    FPdata$bf_profile <- bf_profile
    FPdata$pvalue <- pvalue
    FPdata$tvalue <- tvalue
    FPdata$null_sim <- null_sim
    FPdata$wg_top10_signaltotals <- topten_signal_totals
    FPdata$wg_top10_motifsignal_totals <- topten_motifsignal_totals
    FPdata$wg_signaltotals <- signal_totals
    FPdata$wg_motifsignal_totals <- motifsignal_totals
    ##
    rm(PWM, sites, total_sites, wgsignals, bf_pass_bindingsites, bf_pvalue_pass, bf_profile, pvalue, tvalue, null_sim, topten_signal_totals, topten_motifsignal_totals, signal_totals, motifsignal_totals)
    cat("Saving signals data for gene", gene_name, "motif", a, "at path:", signal_path, "\n")
    save.image(file = signal_path)
    rm(FPdata)
    gc()
  }
}

############

#### Save
cat("Finalizing...\n")
file.create(file = out_filepath)

#### Cleanup workspace
cat("Cleaning workspace...", "\n")
rm(list=ls())
gc()


#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)


#### Library loads
cat("Loading libraries...", "\n")
library(MotifDb)
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Snakemake variables
cat("Loading snakemake variables", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitesfile <- snakemake@input[[3]]
out_filepath <- snakemake@output[[1]]

#### Load workspace
cat("Loading workspace...", "\n")
load(sitesfile)

#### Parameters
cat("Setting input parameters...", "\n")
motif_score <- "95%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens
num_motif <- length(bindingSites)

signalList <- list()

#### Make the footprints
for (a in 1:num_motif){
  
  cat("Generating footprint signal", a, "of ", num_motif,  "\n")
  PWM <- bindingSites[[a]][["PWM"]]
  sites <- bindingSites[[a]][["sites"]]
  sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
  sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
  
  # generate signal
  signal <- factorFootprints(
    # bam file input
    bamfiles = bampath,
    # bai index input
    index = baipath,
    # PWM input
    pfm = PWM,
    # reference genome
    genome = genome,
    # minimum motif matching score
    min.score = motif_score,
    # scope of analysis
    seqlev = scope,
    # bp upstream of motif to analyze
    upstream = upstream,
    # bp downstream of motif to analyze
    downstream = downstream,
    # binding sites
    bindingSites = sites
  )
  
  signalList[[b]] <- signal
  
}

#### Save
cat("Saving at output path: ", out_filepath, "...\n")
save.image(file = out_filepath)

#### Cleanup workspace
cat("Cleaning workspace...", "\n")
rm(list=ls())
gc()


#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("ATACseqQC", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates=TRUE)
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("ChIPpeakAnno", suppressUpdates = TRUE)


#### Library loads
cat("Loading libraries...", "\n")
library(MotifDb)
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)


#### Functions
cat("Creating functions...", "\n")
getMotifPWM <- function(symbol)
{
  
  #### Query MotifDb ####
  # Define string to return only Hsapiens motifs
  organism_rows = grep('Hsapiens', values(MotifDb)$organism, ignore.case = TRUE)
  # Define string for given gene
  gene_symbol_rows = grep(symbol, values(MotifDb)$geneSymbol, ignore.case = TRUE)
  # Get indices for the intersection of gene and organism
  human_gene_rows = intersect(gene_symbol_rows, organism_rows)
  # Pull the PWMs
  # need to make it a list to use unique() function to remove duplicate entries
  gene_motifs <- as.list(MotifDb[human_gene_rows])
  # Get unique motifs
  unique_gene_motifs <- unique(gene_motifs)
  # Get number of motifs
  number_unique_gene_motifs <- length(unique_gene_motifs)
  #### end ####
  
  #### Return PWM ####
  # Make a vector containing all the motifs
  unique_gene_motifs.vector <- c()
  for (a in 1:number_unique_gene_motifs) {
    unique_gene_motifs.vector[a] <- unique_gene_motifs[a]
  }
  
  # return the vector
  return(unique_gene_motifs.vector)
  #### end ####
  
} # end getMotifPWM function


#### Get the PWM for current gene
cat("Fetching PWMs...", "\n")
# pull the current gene symbol from snakemake
gene <- snakemake@wildcards[["gene"]]
# get all human annotated records
PWMlist <- getMotifPWM(gene)


#### Set parameters
cat("Setting input parameters...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
motif_score <- "90%"
upstream <- 100
downstream <- 100
scope <- paste0("chr", c(1:22, "X", "Y"))
genome <- Hsapiens
PWM <- PWMlist[[1]]


#### Generate the footprint
cat("Generating footprint signal...", "\n")
# generate signal
signal <- factorFootprints(
  # bam file input
  bamfiles = bampath,
  # bai index input
  index = baipath,
  # PWM input
  pfm = PWM,
  # reference genome
  genome = genome,
  # minimum motif matching score
  min.score = motif_score,
  # scope of analysis
  seqlev = scope,
  # bp upstream of motif to analyze
  upstream = upstream,
  # bp downstream of motif to analyze
  downstream = downstream
)


#### Save
cat("Saving...", "\n")
outpath <- snakemake@output[[1]]
save(outpath)




