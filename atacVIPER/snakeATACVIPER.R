
##
source("https://bioconductor.org/biocLite.R")
biocLite("bcellViper", suppressUpdates = TRUE)
biocLite("viper", suppressUpdates = TRUE)


# 1.1 Load VIPER library
library(viper)
# 1.2 Load data
data(bcellViper, package="bcellViper")
# 1.3 Get expression matrix
eset <- exprs(dset)

##
# 1.1 Get a list of all targets of every regulator in the network.
tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
# 1.2 Remove expression data of all the genes which arent in the list identified above.
eset <- eset[rownames(eset) %in% unique(tmp),]
# 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ].
tt <- t(scale(t(eset)))

##
# 2.1 Remove targets in the regulon which are not in the rownames of the expression matrix.
regulon <- lapply(regulon, function(x, genes) {
  filtro <- names(x$tfmode) %in% genes
  x$tfmode <- x$tfmode[filtro]
  if (length(x$likelihood) == length(filtro))
    x$likelihood <- x$likelihood[filtro]
  return(x)
}, genes = rownames(eset))
# 2.2 Define minimum regulon size for filtering (default is 20).
minsize <- 20
# 2.3 Remove regulators with a regulon size below the minsize parameter.
regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= minsize]

##
# 1.1 Get a list of all targets of every regulator in the network.
targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
# 1.2 Create the Mode of Regulation matrix from the regulon object.
mor <- sapply(regulon, function(x, genes) {
  return(x$tfmode[match(genes, names(x$tfmode))])
}, genes = targets)
# 1.2 Create the Weights matrix from the regulon object.
wts <- sapply(regulon, function(x, genes) {
  tmp <- x$likelihood[match(genes, names(x$tfmode))]
  tmp[is.na(match(genes, names(x$tfmode)))] <- NA
  return(tmp/max(tmp, na.rm = T))
}, genes = targets)
# 1.3 For each regulator, assign values of 0 to genes which are not listed as its targets.
mor[is.na(mor)] <- 0
wts[is.na(wts)] <- 0
# 1.4 Scale the columns of the Weights matrix to the sum of the weights.
wtss <- scale(wts, center = FALSE, scale = colSums(wts))

##
# 2.1 Calculate the T2 rank matrix from the expression dataset.
t2 <- apply(tt, 2, rank)/(nrow(tt) + 1)
# This line of code is necessary to match the order of genes
# for the subsequent matrix multiplication steps.
pos <- match(targets, rownames(tt))
# 2.2 Transform T2 ranks to Gaussian values.
t2q <- qnorm(filterRowMatrix(t2, pos))
# 2.3 Matrix multiplication.
sum1 <- t(mor * wtss) %*% t2q


##
# 3.1 Calculate the T1 Rank matrix from the T2 Rank matrix.
t1 <- abs(t2 - 0.5) * 2
t1 <- t1 + (1 - max(t1))/2
# 3.2 Transform T1 ranks to Gaussian values.
t1q <- qnorm(filterRowMatrix(t1, pos))
# 3.3 Matrix multiplication.
sum2 <- t((1 - abs(mor)) * wtss) %*% t1q

##
# 4.1 Extract the signs of the Two-tail enrichment scores
ss <- sign(sum1)
ss[ss == 0] <- 1
# 4.2 Combine the Two-tail and One-tail enrichment score matrices.
sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss


##
# 5.1 For each regulator, calculate an index proportional to the likelihood value
# of all its regulatory interactions.
lwt <- sqrt(colSums(wts^2))
# 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above.
nes <- sum3 * lwt


##
# Remove all variables except for the NES results of 5.2
#rm(list=setdiff(ls(), nes))
# Load VIPER library
library(viper)
# Load data
#data(bcellViper, package="bcellViper")
# Get expression matrix
#eset <- exprs(dset)
# Run VIPER
#nes_original <- viper(eset, regulon, verbose=FALSE)


##
#identical(round(nes, digits=5), round(nes_original, digits=5))



simpleVIPER <- function (eset, regulon, minsize=20)
{
  # Data Preprocessing
  # Step 1 - Filter the expression data.
  # 1.1 Get a list of all targets of every regulator in the network.
  tmp <- c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  # 1.2 Remove expression data of all the genes which arent in the list identified above.
  eset <- eset[rownames(eset) %in% unique(tmp),]
  # 1.3 Scale the expression matrix on rows [i.e. row = (row - mean(row)) / sd(row) ].
  tt <- t(scale(t(eset)))
  # Step 2 - Filter the regulon object (i.e. the regulatory network).
  # 2.1 Remove targets in the regulon which are not in the rownames of the expression matrix.
  regulon <- lapply(regulon, function(x, genes) {
    filtro <- names(x$tfmode) %in% genes
    x$tfmode <- x$tfmode[filtro]
    if (length(x$likelihood) == length(filtro))
      x$likelihood <- x$likelihood[filtro]
    return(x)
  }, genes = rownames(eset))
  # 2.2 Define minimum regulon size for filtering (default is 20).
  # The minsize parameter is specified in the function parameters.
  # 2.3 Remove regulators with a regulon size below the minsize parameter.
  regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= minsize]
  # aREA
  nes <- simpleaREA(tt, regulon)
  # Return VIPER matrix
  return(nes)
}



simpleaREA <- function (tt, regulon)
{
  # Step 1 - Create the Mode of Regulation and Weights matrices.
  # 1.1 Get a list of all targets of every regulator in the network.
  targets <- unique(unlist(lapply(regulon, function(x) names(x$tfmode)), use.names = FALSE))
  # 1.2 Create the Mode of Regulation matrix from the regulon object.
  mor <- sapply(regulon, function(x, genes) {
    return(x$tfmode[match(genes, names(x$tfmode))])
  }, genes = targets)
  # 1.2 Create the Weights matrix from the regulon object.
  wts <- sapply(regulon, function(x, genes) {
    tmp <- x$likelihood[match(genes, names(x$tfmode))]
    tmp[is.na(match(genes, names(x$tfmode)))] <- NA
    return(tmp/max(tmp, na.rm = T))
  }, genes = targets)
  # 1.3 For each regulator, assign values of 0 to genes which are not listed as its targets.
  mor[is.na(mor)] <- 0
  wts[is.na(wts)] <- 0
  # 1.4 Scale the columns of the Weights matrix to the sum of the weights.
  wtss <- scale(wts, center = FALSE, scale = colSums(wts))
  # Step 2 - Calculate the two-tail enrichment scores.
  # 2.1 Calculate the T2 rank matrix from the expression dataset.
  t2 <- apply(tt, 2, rank)/(nrow(tt) + 1)
  # This line of code is necessary to match the order of genes
  # for the subsequent matrix multiplication steps.
  pos <- match(targets, rownames(tt))
  # 2.2 Transform T2 ranks to Gaussian values.
  t2q <- qnorm(filterRowMatrix(t2, pos))
  # 2.3 Matrix multiplication.
  sum1 <- t(mor * wtss) %*% t2q
  # Step 3 - Calculate the one-tail score matrix.
  # 3.1 Calculate the T1 Rank matrix from the T2 Rank matrix.
  t1 <- abs(t2 - 0.5) * 2
  t1 <- t1 + (1 - max(t1))/2
  # 3.2 Get qnorm values
  t1q <- qnorm(filterRowMatrix(t1, pos))
  # 3.3 Matrix multiplication.
  sum2 <- t((1 - abs(mor)) * wtss) %*% t1q
  # Step 4 - Calculate the Three-tail enrichment score.
  # 4.1 Extract the signs of the Two-tail enrichment scores
  ss <- sign(sum1)
  ss[ss == 0] <- 1
  # 4.2 Combine the Two-tail and One-tail enrichment score matrices.
  sum3 <- (abs(sum1) + sum2 * (sum2 > 0)) * ss
  # Step 5 - Calculate the Normalized Enrichment Scores.
  # 5.1 For each regulator, calculate an index proportional to the likelihood value
  # of all its regulatory interactions.
  lwt <- sqrt(colSums(wts^2))
  # 5.2 Adjust 3T enrichment scores proportionally to the weights calculated above.
  nes <- sum3 * lwt
}




x <- simpleVIPER(eset, regulon)

y <- msviper(eset, regulon)

plot(y, cex=.7)
