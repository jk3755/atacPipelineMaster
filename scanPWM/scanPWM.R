
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

## Query the database to return all human annotations
mdbHuman <- query (MotifDb, 'hsapiens')
## Get the list of all unique genes
uniqueGenes <- unique(mdbHuman@elementMetadata@listData[["geneSymbol"]])
## Get the total number of unique annotated genes
numGenes <- length(uniqueGenes)
## Initiate list object to store all motif data
motifData <- list()

##
for (gene in uniqueGenes){
  
    ## Get relevant indices
    geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
    ## Number of redundant motifs for current gene
    numMotifs <- length(geneIdx)
       
    ## Get the unique motifs for current gene
    tempMotifs <- list()
    c <- 1
    
    ##
    for (idx in geneIdx){
      tempMotifs[c] <- mdbHuman@listData[idx]
      c <- c+1} # end for (idx in geneIdx)
    
    ##
    uniqueMotifs <- unique(tempMotifs)
    numUniqueMotifs <- length(uniqueMotifs)
    
    ## Initiate sublists for current gene
    tryCatch({
      com <- paste0("motifData$", gene, " <- list()")
      eval(parse(text = com))},
    error=function(cond){
      message(cond)
      return(NA)},
    finally={}) #end tryCatch()
             
    ## Populate motifData list for current gene
    tryCatch({     
      for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
     },
    error=function(cond){
      message(cond)
      return(NA)},
      finally={}) #end tryCatch()

} # end for (gene in uniqueGenes)


> ## Extract relevant data
  > for (gene in uniqueGenes){
    +   
      +   ## Get relevant indices
      +   geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
      +   ## Number of redundant motifs for current gene
        +   numMotifs <- length(geneIdx)
        +   
          +   ## Get the unique motifs for current gene
          +   tempMotifs <- list()
          +   c <- 1
          +   for (idx in geneIdx){
            +     tempMotifs[c] <- mdbHuman@listData[idx]
            +     c <- c+1
            +   } # end for (idx in geneIdx)
          +   ##
            +   uniqueMotifs <- unique(tempMotifs)
            +   numUniqueMotifs <- length(uniqueMotifs)
            +   
              +   ## Initiate sublists for current gene
              +   tryCatch({
                +     com <- paste0("motifData$", gene, " <- list()")
                +     eval(parse(text = com))},
                +   error=function(cond){
                  +     message(cond)
                  +     return(NA)},
                +   finally={}) #end tryCatch()
            +   
              +   ## Populate motifData list for current gene
              +   tryCatch({
                +     
                  +     for (a in 1:numUniqueMotifs){
                    +       com <- paste0("motifData$", gene, "$motif", a, " <- list()")
                    +       eval(parse(text = com))
                    +     } # end for (a in 1:numUniqueMotifs)
                +     },
                +     error=function(cond){
                  +       message(cond)
                  +       return(NA)},
                +     finally={}) #end tryCatch()
            +   
              + } # end for (gene in uniqueGenes)


> ## Extract relevant data
  > for (gene in uniqueGenes){
    +   
      +   ## Get relevant indices
      +   geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
      +   ## Number of redundant motifs for current gene
        +   numMotifs <- length(geneIdx)
        +   
          +   ## Get the unique motifs for current gene
          +   tempMotifs <- list()
          +   c <- 1
          +   for (idx in geneIdx){
            +     tempMotifs[c] <- mdbHuman@listData[idx]
            +     c <- c+1
            +   } # end for (idx in geneIdx)
          +   ##
            +   uniqueMotifs <- unique(tempMotifs)
            +   numUniqueMotifs <- length(uniqueMotifs)
            +   
              +   ## Initiate sublists for current gene
              +   f <- tryCatch({
                +     com <- paste0("motifData$", gene, " <- list()")
                +     eval(parse(text = com))
                +   },
                +   error=function(cond){
                  +     message(cond)
                  +     return(NA)
                  +   },
                +   finally={}) #end tryCatch()
              +   
                +   # ## Populate motifData list for current gene
                +   # for (a in 1:numUniqueMotifs){
                +   #   
                +   #   ##
                +   #   com <- paste0("motifData$", gene, "$motif", a, " <- list()")
                +   #   eval(parse(text = com))
                +     
                +  # } # end for (a in 1:numMotifs)
                + } # end for (gene in uniqueGenes)


> uniqueMotifs <- unique(tempMotifs)
> View(uniqueMotifs)
> View(uniqueMotifs)
> numUniqueMotifs <- length(uniqueMotifs)
> 
              > 