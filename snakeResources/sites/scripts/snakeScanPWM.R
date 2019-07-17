
#### Library loading
cat("Loading libraries for PWM scan", "\n")
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Biostrings))

##
cat("Setting snakemake variables for PWM scan", "\n")
motifDataPath <- snakemake@input[[1]]
outPath <- snakemake@output[[1]]
currentGene <- snakemake@wildcards[["gene"]]

##
cat("Loading motifData object", "\n")
load(motifDataPath)

##
cat("Creating binding sites object for gene", currentGene, "\n")
com <- paste0("motifs <- motifData$", currentGene)
eval(parse(text = com))

## Set parameters
cat("Setting parameters", "\n")
score <- "95%"
numMotifs <- length(motifs)
genome <- Hsapiens
bindingSites <- list()

##
cat("Found", numMotifs, "unique motifs for gene", currentGene, "scanning for sites with", score, "PWM match", "\n")

## Scan the genome for matches to each unique motif
for (a in 1:numMotifs){
  
  cat("Scanning motif", a, "for gene", currentGene, "\n")
  
  cat("Creating temp sites object", "\n")
  tempSites <- list()
  
  cat("Loading PWM", "\n")
  PWM <- motifs[[a]]
  
  cat("Scanning for sites", "\n")
  sites <- matchPWM(PWM, genome, min.score = score)
  
  cat("Transferring data", "\n")
  tempSites$PWM <- PWM
  tempSites$sites <- sites
  bindingSites[[a]] <- tempSites
  
} # end for (a in 1:numMotifs)

## Save the data
cat("Saving binding sites data at path:", outPath, "\n")
save(bindingSites, file = outPath)