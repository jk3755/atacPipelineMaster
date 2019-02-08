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
num_pwm <- length(PWMlist)
cat("Found", num_pwm, " unique motifs for gene", gene, "\n")
genome <- Hsapiens
bindingSites <- list()

#### Scan
for (b in 1:num_pwm){
  
  cat("Scanning motif", b, "for gene", gene, "\n")
  temp <- list()
  PWM <- PWMlist[[b]]
  sites <- matchPWM(PWM, genome, min.score="90%")
  #
  temp$PWM <- PWM
  temp$sites <- sites
  bindingSites[[b]] <- temp
  
} # end for (b in 1:num_pwm)

#### Save
cat("Saving...", "\n")
outpath <- snakemake@output[[1]]
save(bindingSites, file = outpath)

#### Cleanup
cat("Cleaning workspace...", "\n")
rm(list=ls())
gc()

cat("Finished!", "\n")