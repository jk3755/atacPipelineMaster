#### Generate gene uniqueGenes in snakemake group format and output to text file. Only need to run once
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)


##############################################################################################################################################
##############################################################################################################################################
## Raw TF group rule names
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
  for (idx in geneIdx){
    tempMotifs[c] <- mdbHuman@listData[idx]
    c <- c+1} # end for (idx in geneIdx)
  uniqueMotifs <- unique(tempMotifs)
  numUniqueMotifs <- length(uniqueMotifs)
  ## Initiate sublists for current gene
  tryCatch({
    com <- paste0("motifData$", gene, " <- list()")
    eval(parse(text = com))},
    error=function(cond){
      #message(cond)
      return(NA)},
    finally={}) #end tryCatch()
  ## Populate motifData list for current gene
  tryCatch({     
    for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
  },
  error=function(cond){
    #message(cond)
    return(NA)},
  finally={}) #end tryCatch()
  
} # end for (gene in uniqueGenes)
numGenes <- length(motifData)
uniqueGenes <- names(motifData)
strings <- c()
a <- 1 # count for genes
b <- 1 # string index
while (a <= numGenes){

  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4
  h <- a+5
  i <- a+6
  j <- a+7
  k <- a+8
  l <- a+9
  m <- a+10
  n <- a+11
  o <- a+12
  p <- a+13
  q <- a+14
  r <- a+15
  s <- a+16
  t <- a+17
  u <- a+18
  v <- a+19
  
  tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[c], ".rawFPanalysis.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[d], ".rawFPanalysis.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[e], ".rawFPanalysis.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[f], ".rawFPanalysis.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[g], ".rawFPanalysis.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[h], ".rawFPanalysis.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[i], ".rawFPanalysis.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[j], ".rawFPanalysis.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[k], ".rawFPanalysis.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[l], ".rawFPanalysis.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[m], ".rawFPanalysis.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[n], ".rawFPanalysis.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[o], ".rawFPanalysis.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[p], ".rawFPanalysis.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[q], ".rawFPanalysis.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[r], ".rawFPanalysis.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[s], ".rawFPanalysis.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[t], ".rawFPanalysis.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[u], ".rawFPanalysis.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", uniqueGenes[v], ".rawFPanalysis.bamcopy20.done", "' ")
  
  strings[b] <- paste0(
                      "rule rawTF_group",
                      b,
                      ":\n",
                      "\tinput:\n\t\t",
                      tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
                      tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
                      tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
                      tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
                      "output:\n\t\t",
                      "'{path}footprints/operations/groups/{mergedsample}.rawTF.group", b, ".done'\n\t",
                      "shell:\n\t\t",
                      "'touch {output}'"
                      )
  a <- a+20
  b <- b+1
  
}

## Write the file
outPath <- "C:\\Users\\jsk33\\Desktop\\n.txt"
write.table(
            strings,
            file = outPath,
            quote = FALSE,
            sep = ",",
            eol = "\n",
            row.names = FALSE,
            col.names = FALSE)

##############################################################################################################################################
##############################################################################################################################################
## Raw TF run rule names

strings <- c()
for (b in 1:62){
  
  strings[b] <- paste0(
    "rule run_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    "'snu61/wt01/footprints/operations/groups/SNU61-WT-01.rawTF.group", b, ".done'"
  )}

## Write the file
##
outPath <- "C:\\Users\\jsk33\\Desktop\\n2.txt"
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)



##############################################################################################################################################
##############################################################################################################################################
## Parse TF group names

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
      #message(cond)
      return(NA)},
    finally={}) #end tryCatch()
  ## Populate motifData list for current gene
  tryCatch({     
    for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
  },
  error=function(cond){
    #message(cond)
    return(NA)},
  finally={}) #end tryCatch()
} # end for (gene in uniqueGenes)
numGenes <- length(motifData)
uniqueGenes <- names(motifData)
strings <- c()
##
a <- 1 # count for genes
b <- 1 # string index
##
while (a <= numGenes){
  ##
  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4
  h <- a+5
  i <- a+6
  j <- a+7
  k <- a+8
  l <- a+9
  m <- a+10
  n <- a+11
  o <- a+12
  p <- a+13
  q <- a+14
  r <- a+15
  s <- a+16
  t <- a+17
  u <- a+18
  v <- a+19
  ##
  tmp1 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[c], ".parseFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[d], ".parseFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[e], ".parseFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[f], ".parseFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[g], ".parseFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[h], ".parseFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[i], ".parseFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[j], ".parseFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[k], ".parseFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[l], ".parseFP.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[m], ".parseFP.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[n], ".parseFP.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[o], ".parseFP.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[p], ".parseFP.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[q], ".parseFP.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[r], ".parseFP.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[s], ".parseFP.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[t], ".parseFP.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[u], ".parseFP.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/parse/{mergedsample}.", uniqueGenes[v], ".parseFP.bamcopy20.done", "' ")
  ##
  strings[b] <- paste0(
    "rule parseTF_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.parseTF.group", b, ".done'\n\t",
    "shell:\n\t\t",
    "'touch {output}'"
  )
  a <- a+20
  b <- b+1
  
} # end while (a <= numGenes)


## Write the file
outPath <- "C:\\Users\\jsk33\\Desktop\\names.txt"
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)

##############################################################################################################################################
##############################################################################################################################################
## Parse TF run rule names

strings <- c()
for (b in 1:62){
  
  strings[b] <- paste0(
    "rule parse_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    "'snu61/wt01/footprints/operations/groups/SNU61-WT-01.parseTF.group", b, ".done'"
  )}

## Write the file
##
outPath <- "C:\\Users\\jsk33\\Desktop\\runrule.txt"
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)

##############################################################################################################################################
##############################################################################################################################################
## Raw TF group rule names
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
      #message(cond)
      return(NA)},
    finally={}) #end tryCatch()
  ## Populate motifData list for current gene
  tryCatch({     
    for (a in 1:numUniqueMotifs){
      com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
      eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
  },
  error=function(cond){
    #message(cond)
    return(NA)},
  finally={}) #end tryCatch()
} # end for (gene in uniqueGenes)
numGenes <- length(motifData)
uniqueGenes <- names(motifData)
strings <- c()
##
a <- 1 # count for genes
b <- 1 # string index
##
while (a <= numGenes){
  ##
  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4
  h <- a+5
  i <- a+6
  j <- a+7
  k <- a+8
  l <- a+9
  m <- a+10
  n <- a+11
  o <- a+12
  p <- a+13
  q <- a+14
  r <- a+15
  s <- a+16
  t <- a+17
  u <- a+18
  v <- a+19
  ##
  tmp1 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[c], ".parseFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[d], ".parseFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[e], ".parseFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[f], ".parseFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[g], ".parseFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[h], ".parseFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[i], ".parseFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[j], ".parseFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[k], ".parseFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[l], ".parseFP.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[m], ".parseFP.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[n], ".parseFP.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[o], ".parseFP.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[p], ".parseFP.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[q], ".parseFP.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[r], ".parseFP.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[s], ".parseFP.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[t], ".parseFP.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[u], ".parseFP.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/{mergedsample}.", uniqueGenes[v], ".parseFP.bamcopy20.done", "' ")
  ##
  strings[b] <- paste0(
    "rule parseTF_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.processFP.group", b, ".done'\n\t",
    "shell:\n\t\t",
    "'touch {output}'"
  )
  a <- a+20
  b <- b+1
  
} # end while (a <= numGenes)


## Write the file
outPath <- "C:\\Users\\jsk33\\Desktop\\names.txt"
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)


