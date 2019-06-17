
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
  tmp1 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[c], ".graphFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[d], ".graphFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[e], ".graphFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[f], ".graphFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[g], ".graphFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[h], ".graphFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[i], ".graphFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[j], ".graphFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[k], ".graphFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[l], ".graphFP.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[m], ".graphFP.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[n], ".graphFP.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[o], ".graphFP.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[p], ".graphFP.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[q], ".graphFP.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[r], ".graphFP.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[s], ".graphFP.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[t], ".graphFP.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[u], ".graphFP.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[v], ".graphFP.bamcopy20.done", "'")
  
  ##
  strings[b] <- paste0(
    "rule graphFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t",
    tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t",
    tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.graphFP.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+20
  b <- b+1
}

#### Groups 40-55
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
  tmp1 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[c], ".graphFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[d], ".graphFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[e], ".graphFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[f], ".graphFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[g], ".graphFP.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[h], ".graphFP.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[i], ".graphFP.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[j], ".graphFP.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[k], ".graphFP.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[l], ".graphFP.bamcopy10.done", "', ")
  
  ##
  strings[b] <- paste0(
    "rule graphFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
    tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.graphFP.group", b, ".done'\n",
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
  tmp1 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[c], ".graphFP.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[d], ".graphFP.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[e], ".graphFP.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[f], ".graphFP.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}footprints/operations/graphs/{mergedsample}.", orderedNames[g], ".graphFP.bamcopy5.done", "', ")
  
  ##
  strings[b] <- paste0(
    "rule graphFP_group",
    b,
    ":\n",
    "\tinput:\n\t\t",
    tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t",
    "output:\n\t\t",
    "'{path}footprints/operations/groups/{mergedsample}.graphFP.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+5
  b <- b+1
}


#### Write the file ####
outPath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\names\\panTFgraph.snakefile"

##
write.table(
            strings,
            file = outPath,
            quote = FALSE,
            sep = ",",
            eol = "\n",
            row.names = FALSE,
            col.names = FALSE)


