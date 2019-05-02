#### MEMORY taps out after group 43
## Starting at group 40, split the groups into smaller batches
##
#### Generate gene uniqueGenes in snakemake group format and output to text file. Only need to run once
#### Library installs
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
#biocLite("Biostrings", suppressUpdates = TRUE)
#biocLite("MotifDb", suppressUpdates = TRUE)

# #### Library loading
# cat("Loading libraries...", "\n")
# suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
# suppressMessages(library(Biostrings))
# suppressMessages(library(MotifDb))
# 
# ## Query the database to return all human annotations
# mdbHuman <- query (MotifDb, 'hsapiens')
# ## Get the list of all unique genes
# uniqueGenes <- unique(mdbHuman@elementMetadata@listData[["geneSymbol"]])
# ## Get the total number of unique annotated genes
# numGenes <- length(uniqueGenes)
# ## Initiate list object to store all motif data
# motifData <- list()
# 
# ##
# for (gene in uniqueGenes){
#   
#   ## Get relevant indices
#   geneIdx <- which(mdbHuman@elementMetadata@listData[["geneSymbol"]] == gene)
#   ## Number of redundant motifs for current gene
#   numMotifs <- length(geneIdx)
#   
#   ## Get the unique motifs for current gene
#   tempMotifs <- list()
#   c <- 1
#   
#   ##
#   for (idx in geneIdx){
#     tempMotifs[c] <- mdbHuman@listData[idx]
#     c <- c+1} # end for (idx in geneIdx)
#   
#   ##
#   uniqueMotifs <- unique(tempMotifs)
#   numUniqueMotifs <- length(uniqueMotifs)
#   
#   ## Initiate sublists for current gene
#   tryCatch({
#     com <- paste0("motifData$", gene, " <- list()")
#     eval(parse(text = com))},
#     error=function(cond){
#       #message(cond)
#       return(NA)},
#     finally={}) #end tryCatch()
#   
#   ## Populate motifData list for current gene
#   tryCatch({     
#     for (a in 1:numUniqueMotifs){
#       com <- paste0("motifData$", gene, "$motif", a, " <- uniqueMotifs[[a]]")
#       eval(parse(text = com))} # end for (a in 1:numUniqueMotifs)
#   },
#   error=function(cond){
#     #message(cond)
#     return(NA)},
#   finally={}) #end tryCatch()
#   
# } # end for (gene in uniqueGenes)
# 
# 
# uniqueGenes <- names(motifData)





####################################################
namePath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\scanPWM\\bindingSitesSizeOrdered.txt"
#namePath <- "C:\\Users\\Jordan\\Documents\\git\\atacPipelineMaster\\scanPWM\\bindingSitesSizeOrdered.txt"
orderedNames <- readLines(namePath)
numGenes <- length(orderedNames)

##
strings <- c()
##
a <- 1 # count for genes
b <- 1 # string index (group)
##########################################################
# groups 1-40
while (b <= 40){
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
  tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[h], ".rawFPanalysis.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[i], ".rawFPanalysis.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[j], ".rawFPanalysis.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[k], ".rawFPanalysis.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[l], ".rawFPanalysis.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[m], ".rawFPanalysis.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[n], ".rawFPanalysis.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[o], ".rawFPanalysis.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[p], ".rawFPanalysis.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[q], ".rawFPanalysis.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[r], ".rawFPanalysis.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[s], ".rawFPanalysis.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[t], ".rawFPanalysis.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[u], ".rawFPanalysis.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[v], ".rawFPanalysis.bamcopy20.done", "'")
  
  strings[b] <- paste0(
    "rule rawFPanalysis_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'"
  )
  a <- a+20
  b <- b+1
}

############################################################################################
# groups 40-55
while (b <= 55){
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
  
  ##
  tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[h], ".rawFPanalysis.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[i], ".rawFPanalysis.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[j], ".rawFPanalysis.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[k], ".rawFPanalysis.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[l], ".rawFPanalysis.bamcopy10.done", "', ")
  
  strings[b] <- paste0(
    "rule rawFPanalysis_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'"
  )
  a <- a+10
  b <- b+1
}
###########################################################################################
##
while (a <= 1229){
  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4

  ##
  tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
  
  strings[b] <- paste0(
    "rule rawFPanalysis_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'"
  )
  a <- a+5
  b <- b+1
}
###########################################################################################
## Write the file
outPath <- "C:\\Users\\jsk33\\Desktop\\test.txt"
write.table(
            strings,
            file = outPath,
            quote = FALSE,
            sep = ",",
            eol = "\n",
            row.names = FALSE,
            col.names = FALSE)


