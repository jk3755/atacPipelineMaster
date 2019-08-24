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
outpathdone <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
genename <- snakemake@wildcards[["gene"]]
dirpath <- snakemake@wildcards[["path"]]

##
cat("Determining number of motifs...", "\n")
load(sitespath)
num_motifs <- length(bindingSites)
cat("Found ", num_motifs, " motifs", "\n")

for (x in 1:num_motifs){
  
  graphpath <- paste0(dirpath, "footprints/graphs/", samplename, ".", genename, ".", "motif", x, ".svg")
  cat("Output path for signal object: ", graphpath, "\n")
  
  if (file.exists(graphpath) == TRUE){
    cat("File already exists, skipping...", "\n")
    next
    
  } else {
    
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

    load(mergein)
    sigs <- merged_signals[["signal"]]

    ##
    PWM <- bindingSites[[x]][["PWM"]]
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
    
    }}
    

##
cat("Finished...", "\n")
file.create(outpathdone)