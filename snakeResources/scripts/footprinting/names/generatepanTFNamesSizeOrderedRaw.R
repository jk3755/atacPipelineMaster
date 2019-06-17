
##
namePath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\snakeResources\\scripts\\footprinting\\names\\bindingSitesSizeOrdered.txt"
#namePath <- "C:\\Users\\Jordan\\Documents\\git\\atacPipelineMaster\\panTF\\names\\bindingSitesSizeOrdered.txt"

##
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
  tmp1 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
  tmp2 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
  tmp3 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
  tmp4 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
  tmp5 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
  tmp6 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[h], ".rawFPanalysis.bamcopy6.done", "', ")
  tmp7 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[i], ".rawFPanalysis.bamcopy7.done", "', ")
  tmp8 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[j], ".rawFPanalysis.bamcopy8.done", "', ")
  tmp9 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[k], ".rawFPanalysis.bamcopy9.done", "', ")
  tmp10 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[l], ".rawFPanalysis.bamcopy10.done", "', ")
  tmp11 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[m], ".rawFPanalysis.bamcopy11.done", "', ")
  tmp12 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[n], ".rawFPanalysis.bamcopy12.done", "', ")
  tmp13 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[o], ".rawFPanalysis.bamcopy13.done", "', ")
  tmp14 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[p], ".rawFPanalysis.bamcopy14.done", "', ")
  tmp15 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[q], ".rawFPanalysis.bamcopy15.done", "', ")
  tmp16 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[r], ".rawFPanalysis.bamcopy16.done", "', ")
  tmp17 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[s], ".rawFPanalysis.bamcopy17.done", "', ")
  tmp18 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[t], ".rawFPanalysis.bamcopy18.done", "', ")
  tmp19 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[u], ".rawFPanalysis.bamcopy19.done", "', ")
  tmp20 <- paste0("'{path}operations/footprints/raw/{sample}.", orderedNames[v], ".rawFPanalysis.bamcopy20.done", "'")
  
  ##
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
    "'{path}operations/footprints/groups/raw/{sample}.rawFPanalysis.group", b, ".done'\n",
    "\tshell:\n\t\t",
    "'touch {output}'")
  
  ##
  a <- a+20
  b <- b+1
}

# #### Groups 40-45
# while (b <= 45){
#   
#   ##
#   c <- a
#   d <- a+1
#   e <- a+2
#   f <- a+3
#   g <- a+4
#   h <- a+5
#   i <- a+6
#   j <- a+7
#   k <- a+8
#   l <- a+9
#   
#   ##
#   tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
#   tmp2 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[d], ".rawFPanalysis.bamcopy2.done", "', ")
#   tmp3 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[e], ".rawFPanalysis.bamcopy3.done", "', ")
#   tmp4 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[f], ".rawFPanalysis.bamcopy4.done", "', ")
#   tmp5 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[g], ".rawFPanalysis.bamcopy5.done", "', ")
#   tmp6 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[h], ".rawFPanalysis.bamcopy6.done", "', ")
#   tmp7 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[i], ".rawFPanalysis.bamcopy7.done", "', ")
#   tmp8 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[j], ".rawFPanalysis.bamcopy8.done", "', ")
#   tmp9 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[k], ".rawFPanalysis.bamcopy9.done", "', ")
#   tmp10 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[l], ".rawFPanalysis.bamcopy10.done", "', ")
#   
#   ##
#   strings[b] <- paste0(
#     "rule rawFPanalysis_group",
#     b,
#     ":\n",
#     "\tinput:\n\t\t",
#     tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t",
#     tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t",
#     "output:\n\t\t",
#     "'{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group", b, ".done'\n",
#     "\tshell:\n\t\t",
#     "'touch {output}'")
#   
#   ##
#   a <- a+10
#   b <- b+1
# }

# 
# ## Remaining
# while (a <= 1229){
#   
#   ##
#   c <- a
#   ##
#   tmp1 <- paste0("'{path}footprints/operations/raw/{mergedsample}.", orderedNames[c], ".rawFPanalysis.bamcopy1.done", "', ")
# 
#   ##
#   strings[b] <- paste0(
#     "rule rawFPanalysis_group",
#     b,
#     ":\n",
#     "\tinput:\n\t\t",
#     tmp1, "\n\t",
#     "output:\n\t\t",
#     "'{path}footprints/operations/groups/{mergedsample}.rawFPanalysis.group", b, ".done'\n",
#     "\tshell:\n\t\t",
#     "'touch {output}'")
#   
#   ##
#   a <- a+1
#   b <- b+1
# }


#### Write the file ####
#outPath <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\names\\panTFraw.snakefile"
#woutPath <- "C:\\Users\\Jordan\\Desktop\\names.txt"
outPath <- "C:\\Users\\jsk33\\Desktop\\names.txt"

##
write.table(
            strings,
            file = outPath,
            quote = FALSE,
            sep = ",",
            eol = "\n",
            row.names = FALSE,
            col.names = FALSE)


