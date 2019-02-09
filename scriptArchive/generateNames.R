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

write.table(strings, file = "/home/ubuntu2/atac/PWM/scripts/names.txt", quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)



###########


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
  
  tmp1 <- paste0("'h508/H508-WT-01.", names[c], ".done.txt", "', ")
  tmp2 <- paste0("'h508/H508-WT-01.", names[d], ".done.txt", "', ")
  tmp3 <- paste0("'h508/H508-WT-01.", names[e], ".done.txt", "', ")
  tmp4 <- paste0("'h508/H508-WT-01.", names[f], ".done.txt", "', ")
  tmp5 <- paste0("'h508/H508-WT-01.", names[g], ".done.txt", "', ")
  tmp6 <- paste0("'h508/H508-WT-01.", names[h], ".done.txt", "', ")
  tmp7 <- paste0("'h508/H508-WT-01.", names[i], ".done.txt", "', ")
  tmp8 <- paste0("'h508/H508-WT-01.", names[j], ".done.txt", "', ")
  tmp9 <- paste0("'h508/H508-WT-01.", names[k], ".done.txt", "', ")
  tmp10 <- paste0("'h508/H508-WT-01.", names[l], ".done.txt", "', ")
  tmp11 <- paste0("'h508/H508-WT-01.", names[c], ".done.txt", "', ")
  tmp12 <- paste0("'h508/H508-WT-01.", names[d], ".done.txt", "', ")
  tmp13 <- paste0("'h508/H508-WT-01.", names[e], ".done.txt", "', ")
  tmp14 <- paste0("'h508/H508-WT-01.", names[f], ".done.txt", "', ")
  tmp15 <- paste0("'h508/H508-WT-01.", names[g], ".done.txt", "', ")
  tmp16 <- paste0("'h508/H508-WT-01.", names[h], ".done.txt", "', ")
  tmp17 <- paste0("'h508/H508-WT-01.", names[i], ".done.txt", "', ")
  tmp18 <- paste0("'h508/H508-WT-01.", names[j], ".done.txt", "', ")
  tmp19 <- paste0("'h508/H508-WT-01.", names[k], ".done.txt", "', ")
  tmp20 <- paste0("'h508/H508-WT-01.", names[l], ".done.txt", "'")
  
  strings[b] <- paste0("rule group", b, ":\n", "\tinput:\n\t\t", tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19, tmp20)
  
  a <- a+20
  b <- b+1
  
}

write.table(strings, file = "/home/ubuntu2/atac/PWM/scripts/names2.txt", quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)

