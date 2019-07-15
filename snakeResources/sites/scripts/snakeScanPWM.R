
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

#### Library loading
cat("Loading libraries...", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))
suppressMessages(library(MotifDb))

##
cat("Setting snakemake vars...", "\n")
motifDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
rdataPath <- gsub("operations", "data", outPath)
rdataPath <- gsub("PWMscan.done", "bindingSites.Rdata", rdataPath)
currentGene <- snakemake@wildcards[["gene"]]

##
cat("Loading motifData object...", "\n")
load(motifDataPath)

##
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))

## Set parameters
cat("Using default PWM matching score of 90%", "\n")
score <- "90%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()
cat("Found", numMotifs, " unique motifs for gene", currentGene, "scanning for sites with", score, "PWM match", "\n")

## Scan the genome for matches to each unique motif
for (a in 1:numMotifs){
  cat("Scanning motif", a, "for gene", currentGene, "\n")
  tempSites <- list()
  PWM <- motifs[[a]]
  sites <- matchPWM(PWM, genome, min.score=score)
  tempSites$PWM <- PWM
  tempSites$sites <- sites
  bindingSites[[a]] <- tempSites
} # end for (a in 1:numMotifs)

## Save the data
cat("Saving data...", "\n")
save(bindingSites, file = rdataPath)
file.create(outPath)

rm(list=ls())
gc()
cat("Finished scanning!", "\n")