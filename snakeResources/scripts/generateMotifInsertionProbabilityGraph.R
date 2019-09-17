#### Encase all code in a tryCatch block, so if anything unexpected goes wrong, pipeline will still run ####
tryCatch({
  
  #### Set snakemake variables ####
  cat("Setting snakemake variables", "\n")
  insertionMatrixPath <- snakemake@input[[1]]
  functionSourcePath <- snakemake@input[[2]]
  svgOutPath <- snakemake@output[[1]]
  sampleName <- snakemake@wildcards[["sample"]]
  repNum <- snakemake@wildcards[["repnum"]]
  refGenome <- snakemake@wildcards[["refgenome"]]
  geneName <- snakemake@wildcards[["gene"]]
  
  #### Generate plot title ####
  plotTitle <- paste0(sampleName, ".", repNum, ".", refGenome, ".", geneName)
  
  #### Report ####
  cat("Generating motif insertion probability graph", "\n")
  cat("Insertion matrix path:", insertionMatrixPath, "\n")
  cat("Source path for functions:", functionSourcePath, "\n")
  cat("SVG output path:", svgOutPath, "\n")
  cat("Plot title:", plotTitle, "\n")
  
  #### Load libraries ####
  cat("Loading libraries", "\n")
  suppressMessages(library(grid))
  suppressMessages(library(gridSVG))
  
  #### Source functions ####
  cat("Loading functions from path:", functionSourcePath, "\n")
  source(functionSourcePath)
  
  #### Load insertion matrix data ####
  cat("Loading insertion matrix and binding sites", "\n")
  load(insertionMatrixPath)
  
  #### Prepare data ####
  cat("Preparing data", "\n")
  bindingSites <- insertionMatrixData[["bindingSites"]]
  insertionMatrix <- insertionMatrixData[["insertionMatrix"]]
  numBindingSites <- length(bindingSites)
  numMatrixRows <- length(insertionMatrix[,1])
  maxMotifWidth <- max(bindingSites@ranges@width)
  
  #### Generate the plot ####
  cat("Generating and saving plot", "\n")
  plotInsertionSiteProbability(
                              insertionMatrix = insertionMatrix,
                              motifStart = 250,
                              motifWidth = maxMotifWidth,
                              plotTitle = plotTitle,
                              svgOutPath = svgOutPath
                              )

}, finally = {
})

