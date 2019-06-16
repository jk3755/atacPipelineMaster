
####
namePath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\names\\bindingSitesSizeOrdered.txt"
orderedNames <- readLines(namePath)
numGenes <- length(orderedNames)
strings <- c()
a <- 1 # count for genes
b <- 1 # string index (group)

#### Groups 1-40
while (b <= 40){
  
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
  tmp1 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[c], ".processFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[d], ".processFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[e], ".processFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[f], ".processFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[g], ".processFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[h], ".processFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[i], ".processFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[j], ".processFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[k], ".processFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[l], ".processFP.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[m], ".processFP.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[n], ".processFP.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[o], ".processFP.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[p], ".processFP.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[q], ".processFP.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[r], ".processFP.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[s], ".processFP.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[t], ".processFP.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[u], ".processFP.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[v], ".processFP.bamcopy20.done", "'")
  
  ##
  strings[b] <- paste0(
    "rule processFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.processFP.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+20
  b <- b+1
}


## Groups 40-55
while (b <= 55){
  
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
  
  ##
  tmp1 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[c], ".processFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[d], ".processFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[e], ".processFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[f], ".processFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[g], ".processFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[h], ".processFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[i], ".processFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[j], ".processFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[k], ".processFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[l], ".processFP.bamcopy10.done", "', ")
  
  ##
  strings[b] <- paste0(
    "rule processFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.processFP.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+10
  b <- b+1
}

## Remaining
while (a <= 1229){
  
  ##
  c <- a
  d <- a+1
  e <- a+2
  f <- a+3
  g <- a+4

  ##
  tmp1 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[c], ".processFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[d], ".processFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[e], ".processFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[f], ".processFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/processed/{mergedsample}.", orderedNames[g], ".processFP.bamcopy5.done", "', ")
  
  ##
  strings[b] <- paste0(
    "rule processFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.processFP.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+5
  b <- b+1
}

#### Write the file ####
outPath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\names\\panTFnameProcess.snakefile"

##
write.table(
  strings,
  file = outPath,
  quote = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = FALSE)
