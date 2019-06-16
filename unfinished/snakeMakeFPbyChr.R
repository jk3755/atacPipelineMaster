##
cat("Loading libraries...", "\n")
#BiocManager::install("ATACseqQC")
#BiocManager::install("MotifDb")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("Rsamtools")
#BiocManager::install("ChIPpeakAnno")
#BiocManager::install("GenomicAlignments")
#BiocManager::install("BiocGenerics")
#BiocManager::install("parallel")
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
  signalpath <- paste0(sample_path, "footprints/data/bychr/", samplename, ".", genename, ".", "motif", a, ".", current_chr, ".Rdata")
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
  
  cat("Generating signal for ", genename, "motif", a, "chromosome ", current_chr, "\n")
  
  func <- tryCatch({
	sigs <- suppressWarnings(factorFootprints(bamfiles = bampath,
                           index = baipath,
                           bindingSites = sites,
                           pfm = PWM,
                           genome = genome,
                           min.score = motif_score,
                           seqlev = scope,
                           upstream = upstream,
                           downstream = downstream))
  cat("Saving signals data...", "\n")
  save(sigs, file = signalpath)
	},
	error=function(cond){
	message(cond)
	return(NA)
	},
	finally={})
	
  } # end if (file.exists(signalpath))
} # end for (a in 1:num_motif)

#
cat("Finished...", "\n")
file.create(final_outpath)