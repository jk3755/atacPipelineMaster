#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  peakPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  outPath <- snakemake@output[[1]]
  
  #### Report ####
  cat("Generating peak idiogram", "\n")
  cat("Peak file path:", peakPath, "\n")
  cat("Source path for functions:", functionSourcePath, "\n")
  cat("SVG output path:", outPath, "\n")
  
  #### Load libraries ####
  cat("Loading libraries", "\n")
  suppressMessages(library(genomation))
  suppressMessages(library(GenomicRanges))
  
  #### Install/load RIdeogram ####
  ## Cant be loaded with conda/bioconductor
  cat("Installing/Loading RIdeogram package", "\n")
  options(repos = "http://cran.us.r-project.org")
  if (!requireNamespace("RIdeogram", quietly = TRUE))
    install.packages("RIdeogram")
  suppressMessages(library(RIdeogram))
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  ### Load karyotype data ####
  cat("Loading karyotype skeleton from package", "\n")
  data(human_karyotype, package="RIdeogram")
  
  ### Load peaks file ####
  cat("Loading peaks from file:", peakPath, "\n")
  grPeaks <- importBED(peakPath)
  
  ### Convert chr names ####
  ## Need to convert from full "chr1" values to only integer values for this function to work
  cat("Converting chr names to integer values", "\n")
  chrNames <- c(as.character(seqnames(grPeaks)))
  chrNames <- gsub("chr", "", chrNames)
  
  ### Convert peaks to dataframe ####
  cat("Converting peaks to a dataframe", "\n")
  peakData <- data.frame(
    Chr = chrNames,
    Start = grPeaks@ranges@start,
    End = grPeaks@ranges@start + grPeaks@ranges@width,
    Value = grPeaks@elementMetadata@listData[["score"]])
  
  ### Generate idiogram ####
  cat("Generating peaks idiogram", "\n")
  ideogram(
    karyotype = human_karyotype,
    overlaid = peakData,
    output = outPath)
  
}, finally = {
})

