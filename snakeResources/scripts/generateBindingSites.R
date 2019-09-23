#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  peakPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  bindingSitesOutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  sampleRep <- snakemake@wildcards[["repnum"]]
  geneName <- snakemake@wildcards[["gene"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  
  #### Report ####
  cat("Generating binding sites:", "\n")
  cat("Peaks file:", peakPath, "\n")
  cat("Sample name:", sampleName, "\n")
  cat("Sample rep:", sampleRep, "\n")
  cat("Gene name:", geneName, "\n")
  cat("Reference genome used:", refGenome, "\n")
  cat("Filepath for loading functions:", functionSourcePath, "\n")
  cat("Output filepath for binding sites:", bindingSitesOutPath, "\n")
  
  #### Filecheck ####
  cat("Checking if output file already exists at path:", bindingSitesOutPath, "\n")
  if (file.exists(bindingSitesOutPath)){
    cat("Output file already exists. Skipping", "\n")
  } else {
    cat("Output file not found. Processing", "\n")
    
    #### Load libraries ####
    cat("Loading libraries", "\n")
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    suppressMessages(library(Biostrings))
    suppressMessages(library(MotifDb))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(Rsamtools))
    suppressMessages(library(GenomicAlignments))
    suppressMessages(library(ChIPseeker))
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    suppressMessages(library(genomation))
    genome <- Hsapiens
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    #### Source functions ####
    cat("Loading functions from path:", functionSourcePath, "\n")
    source(functionSourcePath)
    
    #### Retrieve the binding sites for the current gene ####
    cat("Retrieving binding sites", "\n")
    tempAllSites <- getAllBindingSites(gene = geneName)
    
    #### Subset the binding sites with the peaks
    cat("Subsetting sites to peak regions", "\n")
    peaksGR <- importBED(peakPath)
    bindingSites <- subsetByOverlaps(tempAllSites, peaksGR)
    numSites <- length(bindingSites@ranges)
    cat("Identified", numSites, "total binding sites", "\n")
    
    #### Save the footprint data ####
    cat("Saving binding sites data", "\n")
    save(bindingSites, file = bindingSitesOutPath)
  }
}, finally = {
})