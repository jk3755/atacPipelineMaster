#source("https://bioconductor.org/biocLite.R")
#biocLite("MotifDb", suppressUpdates = TRUE)
library(MotifDb)
# Define string to return only Hsapiens motifs
organism_rows = grep('Hsapiens', values(MotifDb)$organism, ignore.case = TRUE)
humanPWM <- MotifDb[organism_rows]
names <- unique(humanPWM@elementMetadata@listData[["geneSymbol"]])
num_names <- length(names)
strings <- c()
#

###
outfilePath1 <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\panTFNames.txt"
outfilePath2 <- "C:\\Users\\jsk33\\Documents\\git\\atacPipelineMaster\\panTF\\panTFNames2.txt"
### 


a <- 1 # count for genes
b <- 1 # string index
#
while (a <= 1273){
  
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
  
  tmp1 <- paste0("'", names[c], ".sites.Rdata", "', ")
  tmp2 <- paste0("'", names[d], ".sites.Rdata", "', ")
  tmp3 <- paste0("'", names[e], ".sites.Rdata", "', ")
  tmp4 <- paste0("'", names[f], ".sites.Rdata", "', ")
  tmp5 <- paste0("'", names[g], ".sites.Rdata", "', ")
  tmp6 <- paste0("'", names[h], ".sites.Rdata", "', ")
  tmp7 <- paste0("'", names[i], ".sites.Rdata", "', ")
  tmp8 <- paste0("'", names[j], ".sites.Rdata", "', ")
  tmp9 <- paste0("'", names[k], ".sites.Rdata", "', ")
  tmp10 <- paste0("'", names[l], ".sites.Rdata", "', ")
  tmp11 <- paste0("'", names[c], ".sites.Rdata", "', ")
  tmp12 <- paste0("'", names[d], ".sites.Rdata", "', ")
  tmp13 <- paste0("'", names[e], ".sites.Rdata", "', ")
  tmp14 <- paste0("'", names[f], ".sites.Rdata", "', ")
  tmp15 <- paste0("'", names[g], ".sites.Rdata", "', ")
  tmp16 <- paste0("'", names[h], ".sites.Rdata", "', ")
  tmp17 <- paste0("'", names[i], ".sites.Rdata", "', ")
  tmp18 <- paste0("'", names[j], ".sites.Rdata", "', ")
  tmp19 <- paste0("'", names[k], ".sites.Rdata", "', ")
  tmp20 <- paste0("'", names[l], ".sites.Rdata", "'")
  
  strings[b] <- paste0("rule group", b, ":\n", "\tinput:\n\t\t", tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20)
  
  a <- a+20
  b <- b+1
  
}

write.table(strings, file = outfilePath1, quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)

#################################################


#source("https://bioconductor.org/biocLite.R")
#biocLite("MotifDb", suppressUpdates = TRUE)
library(MotifDb)
# Define string to return only Hsapiens motifs
organism_rows = grep('Hsapiens', values(MotifDb)$organism, ignore.case = TRUE)
humanPWM <- MotifDb[organism_rows]
names <- unique(humanPWM@elementMetadata@listData[["geneSymbol"]])
num_names <- length(names)
strings <- c()
#
a <- 1 # count for genes
b <- 1 # string index
#
while (a <= 1273){
  
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
  
  tmp1 <- paste0("'{path}pantf/operations/{mergedsample}.", names[c], ".parsed.pantf.done.txt", "', ")
  tmp2 <- paste0("'{path}pantf/operations/{mergedsample}.", names[d], ".parsed.pantf.done.txt", "', ")
  tmp3 <- paste0("'{path}pantf/operations/{mergedsample}.", names[e], ".parsed.pantf.done.txt", "', ")
  tmp4 <- paste0("'{path}pantf/operations/{mergedsample}.", names[f], ".parsed.pantf.done.txt", "', ")
  tmp5 <- paste0("'{path}pantf/operations/{mergedsample}.", names[g], ".parsed.pantf.done.txt", "', ")
  tmp6 <- paste0("'{path}pantf/operations/{mergedsample}.", names[h], ".parsed.pantf.done.txt", "', ")
  tmp7 <- paste0("'{path}pantf/operations/{mergedsample}.", names[i], ".parsed.pantf.done.txt", "', ")
  tmp8 <- paste0("'{path}pantf/operations/{mergedsample}.", names[j], ".parsed.pantf.done.txt", "', ")
  tmp9 <- paste0("'{path}pantf/operations/{mergedsample}.", names[k], ".parsed.pantf.done.txt", "', ")
  tmp10 <- paste0("'{path}pantf/operations/{mergedsample}.", names[l], ".parsed.pantf.done.txt", "', ")
  tmp11 <- paste0("'{path}pantf/operations/{mergedsample}.", names[m], ".parsed.pantf.done.txt", "', ")
  tmp12 <- paste0("'{path}pantf/operations/{mergedsample}.", names[n], ".parsed.pantf.done.txt", "', ")
  tmp13 <- paste0("'{path}pantf/operations/{mergedsample}.", names[o], ".parsed.pantf.done.txt", "', ")
  tmp14 <- paste0("'{path}pantf/operations/{mergedsample}.", names[p], ".parsed.pantf.done.txt", "', ")
  tmp15 <- paste0("'{path}pantf/operations/{mergedsample}.", names[q], ".parsed.pantf.done.txt", "', ")
  tmp16 <- paste0("'{path}pantf/operations/{mergedsample}.", names[r], ".parsed.pantf.done.txt", "', ")
  tmp17 <- paste0("'{path}pantf/operations/{mergedsample}.", names[s], ".parsed.pantf.done.txt", "', ")
  tmp18 <- paste0("'{path}pantf/operations/{mergedsample}.", names[t], ".parsed.pantf.done.txt", "', ")
  tmp19 <- paste0("'{path}pantf/operations/{mergedsample}.", names[u], ".parsed.pantf.done.txt", "', ")
  tmp20 <- paste0("'{path}pantf/operations/{mergedsample}.", names[v], ".parsed.pantf.done.txt", "', ")
  
  strings[b] <- paste0("rule PANTF_group", b, ":\n", "\tinput:\n\t\t", tmp1, "\n\t\t", tmp2, "\n\t\t", tmp3, "\n\t\t", tmp4, "\n\t\t", tmp5, "\n\t\t", tmp6, "\n\t\t", tmp7, "\n\t\t", tmp8, "\n\t\t", tmp9, "\n\t\t", tmp10, "\n\t\t", tmp11, "\n\t\t", tmp12, "\n\t\t", tmp13, "\n\t\t", tmp14, "\n\t\t", tmp15, "\n\t\t", tmp16, "\n\t\t", tmp17, "\n\t\t", tmp18, "\n\t\t", tmp19, "\n\t\t", tmp20)
  
  a <- a+20
  b <- b+1
  
}

write.table(strings, file = outfilePath2, quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)


