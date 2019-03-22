##
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(MotifDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))

##
cat("Setting snakemake vars...", "\n")
final_outpath <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
repnum <- snakemake@wildcards[["repnum"]]
reptot <- snakemake@wildcards[["reptot"]]
genename <- snakemake@wildcards[["gene"]]
current_chr <- snakemake@wildcards[["chr"]]
sample_path <- snakemake@wildcards[["path"]]
current_prob <- snakemake@wildcards[["prob"]]

## Set the input files
cat("Setting input file paths...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]

##
cat("Setting parameters...", "\n")
motif_score <- "99%"
upstream <- 100
downstream <- 100
scope <- current_chr
genome <- Hsapiens

##
cat("Loading binding sites and PWM...", "\n")
load(sitespath)
num_motif <- length(bindingSites)

## Only do motif 1
  a <- 1
  #
  signalpath <- paste0(sample_path, "saturation/footprints/data/bychr/", samplename, "-REP", repnum, "of", reptot, ".", current_prob, ".", genename, ".", current_chr, ".Rdata")
  cat("Output path for signal object: ", signalpath, "\n")
  
  #
  if (file.exists(signalpath) == TRUE){
    cat("File already exists, skipping...", "\n")
    next
  } else {
    
    #
    cat("File not found, proceeding...", "\n")
    PWM <- bindingSites[[a]][["PWM"]]
    sites <- bindingSites[[a]][["sites"]]
    #
    
    cat("Pruning sites 1...", "\n")
    sites <- keepStandardChromosomes(sites, pruning.mode="coarse")
    #
    cat("Pruning sites 2...", "\n")
    sites <- keepSeqlevels(sites, scope, pruning.mode="coarse")
    
    ## error handling
    # if current sites object has < 5 sites, the pipeline will crash
    
      #
      cat("Generating signal for ", genename, "motif", a, "chromosome ", current_chr, "\n")
      # generate signal
      sigs <- tryCatch(
        {
          factorFootprints(bamfiles = bampath,
                               index = baipath,
                               bindingSites = sites,
                               pfm = PWM,
                               genome = genome,
                               min.score = motif_score,
                               seqlev = scope,
                               upstream = upstream,
                               downstream = downstream)
      
      # Save the data
      cat("Saving signals data...", "\n")
      save(sigs, file = signalpath)
        },
      error=function(cond){
        message(cond)
        return(NA)
      },
      finally={
      })
  } # end if (file.exists(signalpath))

#
cat("Finished...", "\n")
file.create(final_outpath)