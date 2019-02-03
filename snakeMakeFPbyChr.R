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
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]
final_outpath <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
genename <- snakemake@wildcards[["gene"]]
current_chr <- snakemake@wildcards[["chr"]]
sample_path <- snakemake@wildcards[["path"]]

##
cat("Setting parameters...", "\n")
motif_score <- "95%"
upstream <- 100
downstream <- 100
scope <- current_chr
genome <- Hsapiens

##
cat("Loading binding sites and PWM...", "\n")
load(sitespath)
num_motif <- length(bindingSites)

##
for (a in 1:num_motif){ # An iterator for each unique motif for this gene
  
  #
  signalpath <- paste0(sample_path, "footprints/temp/", samplename, ".", genename, ".", "motif", a, ".", current_chr, ".Rdata")
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
  if ((length(sites@ranges@start)) < 6){
    
    cat("No binding sites detected for current chr. Skipping", "\n")
    next
    
  } else {
  
  #
  cat("Generating signal for ", genename, "motif", a, "chromosome ", current_chr, "\n")
  # generate signal
  sigs <- factorFootprints(bamfiles = bampath,
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
  
  } # end if ((length(sites@ranges@start)) < 2)
  } # end if (file.exists(signalpath))
} # end for (a in 1:num_motif)

#
cat("Finished...", "\n")
file.create(final_outpath)