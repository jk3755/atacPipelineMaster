## Install libraries if necessary
#source("https://bioconductor.org/biocLite.R")
#biocLite("stringi", suppressUpdates = TRUE)
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("ensembldb", suppressUpdates = TRUE)
#biocLite("EnsDb.Hsapiens.v86", suppressUpdates = TRUE)
#biocLite("S4Vectors", suppressUpdates = TRUE)
#biocLite("Repitools", suppressUpdates = TRUE)
#biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
#install.packages("ggplot2")

## Load libraries
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(genomation)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(S4Vectors)
library(Repitools)
library(ggplot2)
library(crayon)
library(Rsamtools)


## Shorten variable name for TxDb database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#### The input files

## First file is provided by Jeremy, contains RNAseq data
inputFile <- "C:\\Users\\jsk33\\Documents\\lab\\atac\\snu61_tfc"
expData <- read.table(inputFile, header=TRUE, sep = ',')

## Second file is called peaks from ATACseq data, in .bed format
bedFile <- "C:\\Users\\jsk33\\Documents\\lab\\atac\\atac\\snu61\\wt01\\peaks\\macs2\\merged\\SNU61-WT-01-merged_global_normalization_peaks.narrowPeak"
## Read in the .bed peaks file
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Keep standard ranges only
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")

## The input bam files, for calculating raw signals
inputBam <- "C:\\Users\\jsk33\\Documents\\lab\\atac\\atac\\snu61\\wt01\\preprocessing\\11repmerged\\SNU61-WT-01-repmerged.bam"


#### Retrieve gene information from TxDb
## Make a gene ID based key for retrieving data for genes of interest from Jeremy data
geneKey <- c(as.character(expData[,3]))
## Use the select method to get mapping between tx_name (UCSC) and gene_id (ENTREZ)
annotData <- select(txdb, keys = geneKey, columns = "TXNAME", keytype = "GENEID")
## Make a vector with the tx_names
txNames <- c(annotData[,2])
## Get promoter regions, -200 bp from TSS (can be adjusted) for all transcripts from TxDb
promoters <- promoters(txdb, upstream = 200, downstream = 0)
## Trim the GRanges object to keep standard entries only
promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
## Subset the promoters GRanges object using the generated index
## Note that with multiple transcript variants, this number will be much higher than the gene_id list
promoterData <- promoters[promoters$tx_name %in% txNames]


#### Each hit may refer to multiple transcript variants. Subset based on the first hit for a given gene only
newGeneID <- c()
newTxName <- promoterData@elementMetadata@listData[["tx_name"]]
for (a in 1:length(newTxName)){
  idx <- which(newTxName[a] == annotData[,2])
  newGeneID[a] <- annotData[idx,1]}
## Check how many unique genes there are now
length(unique(newGeneID))
## Add the GeneIDs to the promoters GRanges
promoterData@elementMetadata@listData[["gene_id"]] <- c(newGeneID)


#### There may be duplicate entries for gene IDs in the promoterData Granges, so lets keep only one for each
unIdx <- c()
for (b in 1:length(geneKey)){
  tempIdx <- c(which(promoterData@elementMetadata@listData[["gene_id"]] == geneKey[b]))
  unIdx[b] <- tempIdx[1]}

## As a sanity check
length(unique(unIdx))
## Remove NAs
unIdx <- unIdx[which(is.na(unIdx) == FALSE)]
## Remove duplicate IDS
unIdx <- unique(unIdx)
## As a sanity check
length(unique(unIdx))


## Subset for unique only
uniquePromoters <- promoterData[unIdx]
## Check number of entries
length(unique(uniquePromoters@elementMetadata@listData[["gene_id"]]))


#### Find which gene promoters are accessible (have a peak) using the ATACseq data
overlaps <- findOverlaps(snu61Peaks, uniquePromoters)
## How many of the overlaps referend unique annotated promoters?
idxuq <- unique(overlaps@to)
## Get a GRanges of the genes that have a peak
genePeaks <- uniquePromoters[idxuq]
## How many of these are unique GeneIDs? 
length(unique(genePeaks@elementMetadata@listData[["gene_id"]]))
## Get the list of ENTREZ IDs for genes in this set with active promoters
activeGenePromoters <- unique(genePeaks@elementMetadata@listData[["gene_id"]])


#### Annotate all unique genes as having a peak ("RED") or not ("BLACK")
idx2 <- which(expData[,3] %in% activeGenePromoters)
expData$peak <- "BLACK"
expData$peak[idx2] <- "RED"
## Change -inf values to -11
idx3 <- which(expData[,6] == "-Inf")
expData[idx3,6] <- -11


#### Sort and plot the data
sortedExpData <- expData[order(expData$SNU61_LARGE_INTESTINE_log2, decreasing = FALSE),]
##
plot(sortedExpData[,6], col = sortedExpData$peak)
##
hist(sortedExpData[,6], breaks = 40)


#### Make bins of the genes based on log2 expression
## Plot a histogram/dotplot overlay of the percentage of genes in that bin with/without peak
binData <- list()
binData$"-11"$data <- sortedExpData[which(sortedExpData[,6] <= -10),]
binData$"-10"$data <- sortedExpData[which(sortedExpData[,6] <= -9 & sortedExpData[,6] > -10),]
binData$"-9"$data <- sortedExpData[which(sortedExpData[,6] <= -8 & sortedExpData[,6] > -9),]
binData$"-8"$data <- sortedExpData[which(sortedExpData[,6] <= -7 & sortedExpData[,6] > -8),]
binData$"-7"$data <- sortedExpData[which(sortedExpData[,6] <= -6 & sortedExpData[,6] > -7),]
binData$"-6"$data <- sortedExpData[which(sortedExpData[,6] <= -5 & sortedExpData[,6] > -6),]
binData$"-5"$data <- sortedExpData[which(sortedExpData[,6] <= -4 & sortedExpData[,6] > -5),]
binData$"-4"$data <- sortedExpData[which(sortedExpData[,6] <= -3 & sortedExpData[,6] > -4),]
binData$"-3"$data <- sortedExpData[which(sortedExpData[,6] <= -2 & sortedExpData[,6] > -3),]
binData$"-2"$data <- sortedExpData[which(sortedExpData[,6] <= -1 & sortedExpData[,6] > -2),]
binData$"-1"$data <- sortedExpData[which(sortedExpData[,6] <= 0 & sortedExpData[,6] > -1),]
binData$"0"$data <- sortedExpData[which(sortedExpData[,6] <= 1 & sortedExpData[,6] > 0),]
binData$"1"$data <- sortedExpData[which(sortedExpData[,6] <= 2 & sortedExpData[,6] > 1),]
binData$"2"$data <- sortedExpData[which(sortedExpData[,6] <= 3 & sortedExpData[,6] > 2),]
binData$"3"$data <- sortedExpData[which(sortedExpData[,6] <= 4 & sortedExpData[,6] > 3),]
binData$"4"$data <- sortedExpData[which(sortedExpData[,6] <= 5 & sortedExpData[,6] > 4),]
binData$"5"$data <- sortedExpData[which(sortedExpData[,6] <= 6 & sortedExpData[,6] > 5),]
binData$"6"$data <- sortedExpData[which(sortedExpData[,6] <= 7 & sortedExpData[,6] > 6),]
binData$"7"$data <- sortedExpData[which(sortedExpData[,6] <= 8 & sortedExpData[,6] > 7),]
binData$"8"$data <- sortedExpData[which(sortedExpData[,6] <= 9 & sortedExpData[,6] > 8),]
binData$"9"$data <- sortedExpData[which(sortedExpData[,6] <= 10 & sortedExpData[,6] > 9),]
## Calculate the ratio
binData$"-11"$ratio <- (length(which(binData[["-11"]][["data"]][["peak"]] == "RED")) / length(binData[["-11"]][["data"]][["X"]]))*100
binData$"-10"$ratio <- (length(which(binData[["-10"]][["data"]][["peak"]] == "RED")) / length(binData[["-10"]][["data"]][["X"]]))*100
binData$"-9"$ratio <- (length(which(binData[["-9"]][["data"]][["peak"]] == "RED")) / length(binData[["-9"]][["data"]][["X"]]))*100
binData$"-8"$ratio <- (length(which(binData[["-8"]][["data"]][["peak"]] == "RED")) / length(binData[["-8"]][["data"]][["X"]]))*100
binData$"-7"$ratio <- (length(which(binData[["-7"]][["data"]][["peak"]] == "RED")) / length(binData[["-7"]][["data"]][["X"]]))*100
binData$"-6"$ratio <- (length(which(binData[["-6"]][["data"]][["peak"]] == "RED")) / length(binData[["-6"]][["data"]][["X"]]))*100
binData$"-5"$ratio <- (length(which(binData[["-5"]][["data"]][["peak"]] == "RED")) / length(binData[["-5"]][["data"]][["X"]]))*100
binData$"-4"$ratio <- (length(which(binData[["-4"]][["data"]][["peak"]] == "RED")) / length(binData[["-4"]][["data"]][["X"]]))*100
binData$"-3"$ratio <- (length(which(binData[["-3"]][["data"]][["peak"]] == "RED")) / length(binData[["-3"]][["data"]][["X"]]))*100
binData$"-2"$ratio <- (length(which(binData[["-2"]][["data"]][["peak"]] == "RED")) / length(binData[["-2"]][["data"]][["X"]]))*100
binData$"-1"$ratio <- (length(which(binData[["-1"]][["data"]][["peak"]] == "RED")) / length(binData[["-1"]][["data"]][["X"]]))*100
binData$"0"$ratio <- (length(which(binData[["0"]][["data"]][["peak"]] == "RED")) / length(binData[["0"]][["data"]][["X"]]))*100
binData$"1"$ratio <- (length(which(binData[["1"]][["data"]][["peak"]] == "RED")) / length(binData[["1"]][["data"]][["X"]]))*100
binData$"2"$ratio <- (length(which(binData[["2"]][["data"]][["peak"]] == "RED")) / length(binData[["2"]][["data"]][["X"]]))*100
binData$"3"$ratio <- (length(which(binData[["3"]][["data"]][["peak"]] == "RED")) / length(binData[["3"]][["data"]][["X"]]))*100
binData$"4"$ratio <- (length(which(binData[["4"]][["data"]][["peak"]] == "RED")) / length(binData[["4"]][["data"]][["X"]]))*100
binData$"5"$ratio <- (length(which(binData[["5"]][["data"]][["peak"]] == "RED")) / length(binData[["5"]][["data"]][["X"]]))*100
binData$"6"$ratio <- (length(which(binData[["6"]][["data"]][["peak"]] == "RED")) / length(binData[["6"]][["data"]][["X"]]))*100
binData$"7"$ratio <- (length(which(binData[["7"]][["data"]][["peak"]] == "RED")) / length(binData[["7"]][["data"]][["X"]]))*100
binData$"8"$ratio <- (length(which(binData[["8"]][["data"]][["peak"]] == "RED")) / length(binData[["8"]][["data"]][["X"]]))*100
binData$"9"$ratio <- (length(which(binData[["9"]][["data"]][["peak"]] == "RED")) / length(binData[["9"]][["data"]][["X"]]))*100
##
ratios <- c(binData$"-11"$ratio, binData$"-10"$ratio, binData$"-9"$ratio, binData$"-8"$ratio, binData$"-7"$ratio, binData$"-6"$ratio,
            binData$"-5"$ratio, binData$"-4"$ratio, binData$"-3"$ratio, binData$"-2"$ratio, binData$"-1"$ratio, binData$"0"$ratio,
            binData$"1"$ratio, binData$"2"$ratio, binData$"3"$ratio, binData$"4"$ratio, binData$"5"$ratio, binData$"6"$ratio,
            binData$"7"$ratio, binData$"8"$ratio, binData$"9"$ratio)

##
plot(ratios)


#### Generate data to plot signal intensity (area under peak) against log2 gene expression
## Create the required annotation file to run Rsubread
uniquePromoters@elementMetadata$id <- uniquePromoters@elementMetadata@listData[["gene_id"]]
annotUP <- createAnnotationFile(uniquePromoters)
## Get raw counts for each range
uniquePromoterCounts <- featureCounts(files = inputBam, nthreads = 20, annot.ext = annotUP)


#### FIX THIS LATER ####
expData <- expData[-1114,]
expData <- expData[-1882,]
expData <- expData[-1909,]
expData <- expData[-2478,]
#### FIX THIS LATER ####


## Pull the required log2exp values
log2Matrix <- matrix(data = NA, nrow = length(uniquePromoters@elementMetadata@listData[["gene_id"]]), ncol = 2)
colnames(log2Matrix) <- c("ID", "log2exp")
##
for (x in 1:length(uniquePromoters@elementMetadata@listData[["gene_id"]])){
  
  genetemp <- uniquePromoters@elementMetadata@listData[["gene_id"]][[x]]
  ##
  log2Matrix[x,1] <- genetemp
  log2Matrix[x,2] <- expData[idxtemp,6]}


## Create a matrix to hold the log2 expression and ATACseq signal of the genes
## Rownames are geneID
plotMatrix <- matrix(data = NA, nrow = length(uniquePromoters@ranges@start), ncol = 2)
colnames(plotMatrix) <- c("log2exp", "accessibility")
##
plotMatrix[,1] <- log2Matrix[,2]
## Check that two matrices are in same order and aligned
for (l in 1:2435){
  if (log2Matrix[l,1] != uniquePromoterCounts[["annotation"]][["GeneID"]][[l]]){cat("MISMATCH")}}

## Add the ATACseq signals
plotMatrix[,2] <- uniquePromoterCounts[["counts"]]

#### Sort and plot the data
sortedplotMatrix <- plotMatrix[order(plotMatrix[,1], decreasing = FALSE),]

##
plot(sortedplotMatrix[,1], sortedplotMatrix[,2])
plot(sortedplotMatrix[,2], sortedplotMatrix[,1])

################ above is code i wrote for Jeremy ##########################################################################################################



# This script references: https://biodatascience.github.io/compbio/bioc/ranges.html

##
#source("https://bioconductor.org/biocLite.R")
#biocLite("stringi", suppressUpdates = TRUE)
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("ensembldb", suppressUpdates = TRUE)
#biocLite("EnsDb.Hsapiens.v86", suppressUpdates = TRUE)
#biocLite("S4Vectors", suppressUpdates = TRUE)
#biocLite("Rsubread", suppressUpdates = TRUE)
#biocLite("Repitools", suppressUpdates = TRUE)
#install.packages("ggplot2")
#devtools::install_github("r-lib/crayon")

##
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(S4Vectors)
library(Rsubread)
library(Repitools)
library(ggplot2)
library(crayon)
library(Rsamtools)


# From unmodified genes list, normalize data and compare

# specify input files
input_files <- list()
input_files$snu61_bam1 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S1_S6.lanemerge.dp.bam"
input_files$snu61_bam2 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S2_S4.lanemerge.dp.bam"
input_files$snu61_bam3 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S3_S1.lanemerge.dp.bam"
input_files$ls1034_bam1 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S2_S3.lanemerge.dp.bam"
input_files$ls1034_bam2 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S1_S5.lanemerge.dp.bam"
input_files$ls1034_bam3 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S3_S2.lanemerge.dp.bam"
input_files$h508_bam1 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-2_S2.lanemerge.dp.bam"
input_files$h508_bam2 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-3_S1.lanemerge.dp.bam"
input_files$h508_bam3 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-1_S3.lanemerge.dp.bam"

# find the total reads in all bam files
total_reads <- list()
total_reads$snu61_r1 <- countBam(input_files$snu61_bam1)
total_reads$snu61_r2 <- countBam(input_files$snu61_bam2)
total_reads$snu61_r3 <- countBam(input_files$snu61_bam3)
total_reads$ls1034_r1 <- countBam(input_files$ls1034_bam1)
total_reads$ls1034_r2 <- countBam(input_files$ls1034_bam2)
total_reads$ls1034_r3 <- countBam(input_files$ls1034_bam3)
total_reads$h508_r1 <- countBam(input_files$h508_bam1)
total_reads$h508_r2 <- countBam(input_files$h508_bam2)
total_reads$h508_r3 <- countBam(input_files$h508_bam3)

# calculate normalization factors
replicate_normalization <- matrix(data=NA, ncol=3, nrow=9)
colnames(replicate_normalization) <- c("total reads", "normalization factor", "normalized total reads")
rownames(replicate_normalization) <- c("snu61_r1","snu61_r2","snu61_r3","ls1034_r1","ls1034_r2","ls1034_r3","h508_r1","h508_r2","h508_r3")
#
replicate_normalization[1,1] <- total_reads[["snu61_r1"]][["records"]]
replicate_normalization[2,1] <- total_reads[["snu61_r2"]][["records"]]
replicate_normalization[3,1] <- total_reads[["snu61_r3"]][["records"]]
replicate_normalization[4,1] <- total_reads[["ls1034_r1"]][["records"]]
replicate_normalization[5,1] <- total_reads[["ls1034_r2"]][["records"]]
replicate_normalization[6,1] <- total_reads[["ls1034_r3"]][["records"]]
replicate_normalization[7,1] <- total_reads[["h508_r1"]][["records"]]
replicate_normalization[8,1] <- total_reads[["h508_r2"]][["records"]]
replicate_normalization[9,1] <- total_reads[["h508_r3"]][["records"]]
#
minimums <- list()
minimums$snu61_min <- min(replicate_normalization[1,1],replicate_normalization[2,1],replicate_normalization[3,1]) 
minimums$ls1034_min <- min(replicate_normalization[4,1],replicate_normalization[5,1],replicate_normalization[6,1])
minimums$h508_min <- min(replicate_normalization[7,1],replicate_normalization[8,1],replicate_normalization[9,1])
#
replicate_normalization[1,2] <- replicate_normalization[1,1] / minimums$snu61_min
replicate_normalization[2,2] <- replicate_normalization[2,1] / minimums$snu61_min
replicate_normalization[3,2] <- replicate_normalization[3,1] / minimums$snu61_min
replicate_normalization[4,2] <- replicate_normalization[4,1] / minimums$ls1034_min
replicate_normalization[5,2] <- replicate_normalization[5,1] / minimums$ls1034_min
replicate_normalization[6,2] <- replicate_normalization[6,1] / minimums$ls1034_min
replicate_normalization[7,2] <- replicate_normalization[7,1] / minimums$h508_min
replicate_normalization[8,2] <- replicate_normalization[8,1] / minimums$h508_min
replicate_normalization[9,2] <- replicate_normalization[9,1] / minimums$h508_min
#
replicate_normalization[1,3] <- replicate_normalization[1,1] / replicate_normalization[1,2]
replicate_normalization[2,3] <- replicate_normalization[2,1] / replicate_normalization[2,2]
replicate_normalization[3,3] <- replicate_normalization[3,1] / replicate_normalization[3,2]
replicate_normalization[4,3] <- replicate_normalization[4,1] / replicate_normalization[4,2]
replicate_normalization[5,3] <- replicate_normalization[5,1] / replicate_normalization[5,2]
replicate_normalization[6,3] <- replicate_normalization[6,1] / replicate_normalization[6,2]
replicate_normalization[7,3] <- replicate_normalization[7,1] / replicate_normalization[7,2]
replicate_normalization[8,3] <- replicate_normalization[8,1] / replicate_normalization[8,2]
replicate_normalization[9,3] <- replicate_normalization[9,1] / replicate_normalization[9,2]
#
edb <- EnsDb.Hsapiens.v86
norm_genes <- genes(edb)
rm(edb)
norm_genes <- keepStandardChromosomes(norm_genes, pruning.mode="coarse") # subset to only the standard chromosomes
norm_genes@elementMetadata$id <- norm_genes@elementMetadata@listData[["symbol"]]
annot_norm_gene <- createAnnotationFile(norm_genes)
# get the raw counts for each sample / replicate
gene_count_list <- list()
gene_count_list$snu61_r1 <- featureCounts(files = input_files$snu61_bam1, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$snu61_r2 <- featureCounts(files = input_files$snu61_bam2, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$snu61_r3 <- featureCounts(files = input_files$snu61_bam3, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$ls1034_r1 <- featureCounts(files = input_files$ls1034_bam1, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$ls1034_r2 <- featureCounts(files = input_files$ls1034_bam2, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$ls1034_r3 <- featureCounts(files = input_files$ls1034_bam3, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$h508_r1 <- featureCounts(files = input_files$h508_bam1, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$h508_r2 <- featureCounts(files = input_files$h508_bam2, nthreads = 20, annot.ext = annot_norm_gene)
gene_count_list$h508_r3 <- featureCounts(files = input_files$h508_bam3, nthreads = 20, annot.ext = annot_norm_gene)
#
num_genes <- length(gene_count_list[["snu61_r1"]][["counts"]])
gene_names <- c(gene_count_list[["snu61_r1"]][["annotation"]][["GeneID"]])
#
normalized_counts <- list()
normalized_counts$snu61 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(normalized_counts$snu61) <- c("rep1","rep2","rep3","avg","total")
rownames(normalized_counts$snu61) <- c(gene_names)
normalized_counts$ls1034 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(normalized_counts$ls1034) <- c("rep1","rep2","rep3","avg","total")
rownames(normalized_counts$ls1034) <- c(gene_names)
normalized_counts$h508 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(normalized_counts$h508) <- c("rep1","rep2","rep3","avg","total")
rownames(normalized_counts$h508) <- c(gene_names)
#
for (a in 1:num_genes){
  #
  normalized_counts$snu61[a,1] <- as.integer(gene_count_list[["snu61_r1"]][["counts"]][[a]]/replicate_normalization[1,2])
  normalized_counts$snu61[a,2] <- as.integer(gene_count_list[["snu61_r2"]][["counts"]][[a]]/replicate_normalization[2,2])
  normalized_counts$snu61[a,3] <- as.integer(gene_count_list[["snu61_r3"]][["counts"]][[a]]/replicate_normalization[3,2])
  normalized_counts$snu61[a,4] <- ((normalized_counts$snu61[a,1]+normalized_counts$snu61[a,2]+normalized_counts$snu61[a,3])/3)
  normalized_counts$snu61[a,5] <- normalized_counts$snu61[a,1]+normalized_counts$snu61[a,2]+normalized_counts$snu61[a,3]
  #
  normalized_counts$ls1034[a,1] <- as.integer(gene_count_list[["ls1034_r1"]][["counts"]][[a]]/replicate_normalization[4,2])
  normalized_counts$ls1034[a,2] <- as.integer(gene_count_list[["ls1034_r2"]][["counts"]][[a]]/replicate_normalization[5,2])
  normalized_counts$ls1034[a,3] <- as.integer(gene_count_list[["ls1034_r3"]][["counts"]][[a]]/replicate_normalization[6,2])
  normalized_counts$ls1034[a,4] <- ((normalized_counts$ls1034[a,1]+normalized_counts$ls1034[a,2]+normalized_counts$ls1034[a,3])/3)
  normalized_counts$ls1034[a,5] <- normalized_counts$ls1034[a,1]+normalized_counts$ls1034[a,2]+normalized_counts$ls1034[a,3]
  #
  normalized_counts$h508[a,1] <- as.integer(gene_count_list[["h508_r1"]][["counts"]][[a]]/replicate_normalization[7,2])
  normalized_counts$h508[a,2] <- as.integer(gene_count_list[["h508_r2"]][["counts"]][[a]]/replicate_normalization[8,2])
  normalized_counts$h508[a,3] <- as.integer(gene_count_list[["h508_r3"]][["counts"]][[a]]/replicate_normalization[9,2])
  normalized_counts$h508[a,4] <- ((normalized_counts$h508[a,1]+normalized_counts$h508[a,2]+normalized_counts$h508[a,3])/3)
  normalized_counts$h508[a,5] <- normalized_counts$h508[a,1]+normalized_counts$h508[a,2]+normalized_counts$h508[a,3]}
#
sample_normalization <- matrix(data=NA, ncol=3, nrow=3)
colnames(sample_normalization) <- c("total reads", "normalization factor", "normalized total reads")
rownames(sample_normalization) <- c("snu61", "ls1034", "h508")
#
sample_normalization[1,1] <- replicate_normalization[1,3]+replicate_normalization[2,3]+replicate_normalization[3,3]
sample_normalization[2,1] <- replicate_normalization[4,3]+replicate_normalization[5,3]+replicate_normalization[6,3]
sample_normalization[3,1] <- replicate_normalization[7,3]+replicate_normalization[8,3]+replicate_normalization[9,3]
#
minimums$sample_min <- min(sample_normalization[1,1],sample_normalization[2,1],sample_normalization[3,1])
#
sample_normalization[1,2] <- sample_normalization[1,1] / minimums$sample_min
sample_normalization[2,2] <- sample_normalization[2,1] / minimums$sample_min
sample_normalization[3,2] <- sample_normalization[3,1] / minimums$sample_min
#
sample_normalization[1,3] <- sample_normalization[1,1] /sample_normalization[1,2]
sample_normalization[2,3] <- sample_normalization[2,1] /sample_normalization[2,2]
sample_normalization[3,3] <- sample_normalization[3,1] /sample_normalization[3,2]
#
signal_matrix <- matrix(data=NA,nrow=num_genes,ncol=5)
rownames(signal_matrix) <- c(gene_names)
colnames(signal_matrix) <- c("snu61","ls1034","h508","var","sd")
#
for (b in 1:num_genes){
  #
  signal_matrix[b,1] <- as.integer((normalized_counts[["snu61"]][b,4]/sample_normalization[1,2]))
  signal_matrix[b,2] <- as.integer((normalized_counts[["ls1034"]][b,4]/sample_normalization[2,2]))
  signal_matrix[b,3] <- as.integer((normalized_counts[["h508"]][b,4]/sample_normalization[3,2]))
  signal_matrix[b,4] <- stats::var(signal_matrix[b,1:3])
  signal_matrix[b,5] <- stats::sd(signal_matrix[b,1:3])
}
#
coad_mr <- c("TCF7","MNX1","POU5F1B","ESRRA","CDX2","HNF4A","GMEB2","HOXA3","OVOL1","ASCL2","ZSWIM1","CBFA2T2","TARBP1","KLHL31","GTF2IRD1","ZNF696","NUFIP1","SUPT20H","KAT8","KAT2A","TAF4","RBM39","ZMYND8","L3MBTL1","PLAGL2","NCOA5","ZSWIM3","ADNP","ASXL1")
ind <- which(rownames(signal_matrix) %in% coad_mr)
mr_counts <- signal_matrix[ind,]
#



# Visualizations
# Decide on colors for cell lines here
# Violin plots

# Violin plots
# prep data
violin_data <- data.frame("Signal" = c(mr_counts[,1],mr_counts[,2],mr_counts[,3]),
                          "Sample" = c(rep("snu61", times=29), rep("ls1034",times=29),rep("h508",times=29)))
#
violin_plot <- ggplot(violin_data, aes(x=Sample, y=Signal, fill=Sample)) + 
  geom_violin()
violin_plot


# Change violin plot colors by groups
p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=dose)) +
  geom_violin(trim=FALSE)
p



# Violin for all genes


# Violin plots
# prep data
violin_data_all <- data.frame("Signal" = c(signal_matrix[,1],signal_matrix[,2],signal_matrix[,3]),
                              "Sample" = c(rep("snu61", times=num_genes), rep("ls1034",times=num_genes),rep("h508",times=num_genes)))
#
violin_plot <- ggplot(violin_data_all, aes(x=Sample, y=Signal, fill=Sample)) + 
  geom_violin()
violin_plot


# Heatmaps of the master regulators

mr_names <- rownames(mr_counts)
heatmap_data <- data.frame("Signal" = c(mr_counts[,1],mr_counts[,2],mr_counts[,3]),
                           "Sample" = c(rep("snu61", times=29), rep("ls1034",times=29),rep("h508",times=29)),
                           "Gene" = rep(mr_names))

heatmap.plot <- ggplot(data = heatmap_data, aes(x = Sample, y = Gene)) +
  geom_tile(aes(fill = Signal)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  theme(axis.text.y = element_text(size = 6))

# Preview the heatmap
print(heatmap.plot)


# Heatmaps for individual replicates


#
mr_all_reps <- matrix(data=NA, nrow=29, ncol=9)
rownames(mr_all_reps) <- mr_names
colnames(mr_all_reps) <- c("snu61_r1","snu61_r2","snu61_r3","ls1034_r1","ls1034_r2","ls1034_r3","h508_r1","h508_r2","h508_r3")
#
mr_all_reps[,1] <- normalized_counts[["snu61"]][ind,1]
mr_all_reps[,2] <- normalized_counts[["snu61"]][ind,2]
mr_all_reps[,3] <- normalized_counts[["snu61"]][ind,3]
mr_all_reps[,4] <- normalized_counts[["ls1034"]][ind,1]
mr_all_reps[,5] <- normalized_counts[["ls1034"]][ind,2]
mr_all_reps[,6] <- normalized_counts[["ls1034"]][ind,3]
mr_all_reps[,7] <- normalized_counts[["h508"]][ind,1]
mr_all_reps[,8] <- normalized_counts[["h508"]][ind,2]
mr_all_reps[,9] <- normalized_counts[["h508"]][ind,3]
#
heatmap_data_reps <- data.frame(
  "Signal" = c(mr_all_reps[,1],mr_all_reps[,2],mr_all_reps[,3],mr_all_reps[,4],mr_all_reps[,5],mr_all_reps[,6],mr_all_reps[,7],mr_all_reps[,8],mr_all_reps[,9]),
  "Sample" = c(rep("snu61_r1", times=29), rep("snu61_r2", times=29), rep("snu61_r3", times=29), rep("ls1034_r1",times=29), rep("ls1034_r2",times=29), rep("ls1034_r3",times=29), rep("h508_r1",times=29), rep("h508_r2",times=29), rep("h508_r3",times=29)),
  "Gene" = rep(mr_names))
#
heatmap.plot.reps <- ggplot(data = heatmap_data_reps, aes(x = Sample, y = Gene)) +
  geom_tile(aes(fill = Signal)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  theme(axis.text.y = element_text(size = 6))

# Preview the heatmap
print(heatmap.plot.reps)



#### Copy gene info from ensembl database and convert to Granges

edb <- EnsDb.Hsapiens.v86
genes <- genes(edb)
rm(edb)

# subset to only the standard chromosomes
genes <- keepStandardChromosomes(genes, pruning.mode="coarse")

## Create the modified subsets of annotated genes

#### initialize the modified granges objects ####
granges_list <- list()
genes_200up <- genes
genes_2000up <- genes
genes_200up_gene <- genes
genes_2000up_gene <- genes

#### make the modified granges objects ####
for (a in 1:length(genes)){
  
  ### change positive strand entries
  if (genes[a]@strand@values == "+"){
    
    ### 200 bp upstream only
    genes_200up[a]@ranges@start <- as.integer(genes_200up[a]@ranges@start - 200)
    genes_200up[a]@ranges@width <- as.integer(200)
    ### 2000 bp upstream only
    genes_2000up[a]@ranges@start <- as.integer(genes_2000up[a]@ranges@start - 2000)
    genes_2000up[a]@ranges@width <- as.integer(2000)
    ### 200 bp upstream and gene
    genes_200up_gene[a]@ranges@start <- as.integer(genes_200up_gene[a]@ranges@start - 200)
    genes_200up_gene[a]@ranges@width <- as.integer(genes_200up_gene[a]@ranges@width + 200)
    ### 2000 bp upstream and gene
    genes_2000up_gene[a]@ranges@start <- as.integer(genes_2000up_gene[a]@ranges@start - 2000)
    genes_2000up_gene[a]@ranges@width <- as.integer(genes_2000up_gene[a]@ranges@width + 2000)}
  
  ### change negative strand entries
  if (genes[a]@strand@values == "-"){
    
    ### 200 bp upstream only
    genes_200up[a]@ranges@start <- as.integer(genes_200up[a]@ranges@start + genes_200up[a]@ranges@width)
    genes_200up[a]@ranges@width <- as.integer(200)
    ### 2000 bp upstream only
    genes_2000up[a]@ranges@start <- as.integer(genes_2000up[a]@ranges@start + genes_2000up[a]@ranges@width)
    genes_2000up[a]@ranges@width <- as.integer(2000)
    ### 200 bp upstream and gene
    genes_200up_gene[a]@ranges@width <- as.integer(genes_200up_gene[a]@ranges@width + 200)
    ### 2000 bp upstream and gene
    genes_2000up_gene[a]@ranges@width <- as.integer(genes_2000up_gene[a]@ranges@width + 2000)
  }
}

granges_list$genes_200up <- genes_200up
granges_list$genes_2000up <- genes_2000up
granges_list$genes_200up_gene <- genes_200up_gene
granges_list$genes_2000up_gene <- genes_2000up_gene
rm(genes_200up, genes_2000up, genes_200up_gene, genes_2000up_gene)


#### Positive strand ####
# Get the first instance of a positive and negative strand gene from the original list
firstPlus <- which(genes@strand@values == "+")[1]
cat(black$bold(bgWhite("First plus strand gene found at index:",                (blue$bold(firstPlus)),"\n")))
firstName <- genes@elementMetadata@listData[["gene_name"]][firstPlus]
cat(black$bold(bgWhite("This genes name is:",(blue$bold(firstName)),"\n")))
firstStart <- genes@ranges@start[firstPlus]
cat(black$bold(bgWhite("It starts at position:",(blue$bold(firstStart)),"\n")))
firstWidth <- genes@ranges@width[firstPlus]
cat(black$bold(bgWhite("It has a width of:",(blue$bold(firstWidth)),"\n")))

# 200 bp upstream
name200 <- genes_200up@elementMetadata@listData[["gene_name"]][firstPlus]
start200 <- genes_200up@ranges@start[firstPlus]
width200 <- genes_200up@ranges@width[firstPlus]
cat(red$bold(bgWhite("Checking against:",(blue$bold("200 bp upstream")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name200)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start200)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width200)),"\n")))
if(firstName == name200){sanityChecks[1] = TRUE}else{sanityChecks[1] = FALSE}
if(firstStart == (start200+200)){sanityChecks[2] = TRUE}else{sanityChecks[2] = FALSE}
if(width200 == 200){sanityChecks[3] = TRUE}else{sanityChecks[3] = FALSE}

# 2000 bp upstream
name2000 <- genes_2000up@elementMetadata@listData[["gene_name"]][firstPlus]
start2000 <- genes_2000up@ranges@start[firstPlus]
width2000 <- genes_2000up@ranges@width[firstPlus]
cat(red$bold(bgWhite("Checking against:",(blue$bold("2000 bp upstream")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name2000)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start2000)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width2000)),"\n")))
if(firstName == name2000){sanityChecks[4] = TRUE}else{sanityChecks[4] = FALSE}
if(firstStart == (start2000+2000)){sanityChecks[5] = TRUE}else{sanityChecks[5] = FALSE}
if(width2000 == 2000){sanityChecks[6] = TRUE}else{sanityChecks[6] = FALSE}

# 200 bp upstream and gene
name200gene <- genes_200up_gene@elementMetadata@listData[["gene_name"]][firstPlus]
start200gene <- genes_200up_gene@ranges@start[firstPlus]
width200gene <- genes_200up_gene@ranges@width[firstPlus]
cat(red$bold(bgWhite("Checking against:",(blue$bold("200 bp upstream and gene")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name200gene)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start200gene)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width200gene)),"\n")))
if(firstName == name200gene){sanityChecks[7] = TRUE}else{sanityChecks[7] = FALSE}
if(firstStart == (start200gene+200)){sanityChecks[8] = TRUE}else{sanityChecks[8] = FALSE}
if(width200gene == (200+firstWidth)){sanityChecks[9] = TRUE}else{sanityChecks[9] = FALSE}

# 2000 bp upstream and gene
name2000gene <- genes_2000up_gene@elementMetadata@listData[["gene_name"]][firstPlus]
start2000gene <- genes_2000up_gene@ranges@start[firstPlus]
width2000gene <- genes_2000up_gene@ranges@width[firstPlus]
cat(red$bold(bgWhite("Checking against:",(blue$bold("2000 bp upstream and gene")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name2000gene)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start2000gene)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width2000gene)),"\n")))
if(firstName == name2000gene){sanityChecks[10] = TRUE}else{sanityChecks[10] = FALSE}
if(firstStart == (start2000gene+2000)){sanityChecks[11] = TRUE}else{sanityChecks[11] = FALSE}
if(width2000gene == (2000+firstWidth)){sanityChecks[12] = TRUE}else{sanityChecks[12] = FALSE}

#### Negative strand ####
firstNeg <- which(genes@strand@values == "-")[1]
cat(black$bold(bgWhite("First minus strand gene found at index:",                (blue$bold(firstNeg)),"\n")))
firstName <- genes@elementMetadata@listData[["gene_name"]][firstNeg]
cat(black$bold(bgWhite("This genes name is:",(blue$bold(firstName)),"\n")))
firstStart <- genes@ranges@start[firstNeg]
cat(black$bold(bgWhite("It starts at position:",(blue$bold(firstStart)),"\n")))
firstWidth <- genes@ranges@width[firstNeg]
cat(black$bold(bgWhite("It has a width of:",(blue$bold(firstWidth)),"\n")))

# 200 bp upstream
name200 <- genes_200up@elementMetadata@listData[["gene_name"]][firstNeg]
start200 <- genes_200up@ranges@start[firstNeg]
width200 <- genes_200up@ranges@width[firstNeg]
cat(red$bold(bgWhite("Checking against:",(blue$bold("200 bp upstream")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name200)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start200)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width200)),"\n")))
if(firstName == name200){sanityChecks[13] = TRUE}else{sanityChecks[13] = FALSE}
if(firstStart == (start200-firstWidth)){sanityChecks[14] = TRUE}else{sanityChecks[14] = FALSE}
if(width200 == 200){sanityChecks[15] = TRUE}else{sanityChecks[15] = FALSE}

# 2000 bp upstream
name2000 <- genes_2000up@elementMetadata@listData[["gene_name"]][firstNeg]
start2000 <- genes_2000up@ranges@start[firstNeg]
width2000 <- genes_2000up@ranges@width[firstNeg]
cat(red$bold(bgWhite("Checking against:",(blue$bold("2000 bp upstream")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name2000)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start2000)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width2000)),"\n")))
if(firstName == name2000){sanityChecks[16] = TRUE}else{sanityChecks[16] = FALSE}
if(firstStart == (start2000-firstWidth)){sanityChecks[17] = TRUE}else{sanityChecks[17] = FALSE}
if(width2000 == 2000){sanityChecks[18] = TRUE}else{sanityChecks[18] = FALSE}

# 200 bp upstream and gene
name200gene <- genes_200up_gene@elementMetadata@listData[["gene_name"]][firstNeg]
start200gene <- genes_200up_gene@ranges@start[firstNeg]
width200gene <- genes_200up_gene@ranges@width[firstNeg]
cat(red$bold(bgWhite("Checking against:",(blue$bold("200 bp upstream and gene")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name200gene)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start200gene)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width200gene)),"\n")))
if(firstName == name200gene){sanityChecks[19] = TRUE}else{sanityChecks[19] = FALSE}
if(firstStart == start200gene){sanityChecks[20] = TRUE}else{sanityChecks[20] = FALSE}
if(width200gene == (200+firstWidth)){sanityChecks[21] = TRUE}else{sanityChecks[21] = FALSE}

# 2000 bp upstream and gene
name2000gene <- genes_2000up_gene@elementMetadata@listData[["gene_name"]][firstNeg]
start2000gene <- genes_2000up_gene@ranges@start[firstNeg]
width2000gene <- genes_2000up_gene@ranges@width[firstNeg]
cat(red$bold(bgWhite("Checking against:",(blue$bold("2000 bp upstream and gene")),"\n")))
cat(black$bold(bgWhite("This genes name is:",(blue$bold(name2000gene)),"\n")))
cat(black$bold(bgWhite("It starts at position:",(blue$bold(start2000gene)),"\n")))
cat(black$bold(bgWhite("It has a width of:",(blue$bold(width2000gene)),"\n")))
if(firstName == name2000gene){sanityChecks[22] = TRUE}else{sanityChecks[22] = FALSE}
if(firstStart == start2000gene){sanityChecks[23] = TRUE}else{sanityChecks[23] = FALSE}
if(width2000gene == (2000+firstWidth)){sanityChecks[24] = TRUE}else{sanityChecks[24] = FALSE}


# Trim out of bounds ranges


idx <- GenomicRanges:::get_out_of_bound_index(genes)
if (length(idx) != 0L)
  genes <- genes[-idx]

idx <- GenomicRanges:::get_out_of_bound_index(genes_200up)
if (length(idx) != 0L)
  genes_200up <- genes_200up[-idx]

idx <- GenomicRanges:::get_out_of_bound_index(genes_2000up)
if (length(idx) != 0L)
  genes_2000up <- genes_2000up[-idx]

idx <- GenomicRanges:::get_out_of_bound_index(genes_200up_gene)
if (length(idx) != 0L)
  genes_200up_gene <- genes_200up_gene[-idx]

idx <- GenomicRanges:::get_out_of_bound_index(genes_2000up_gene)
if (length(idx) != 0L)
  genes_2000up_gene <- genes_2000up_gene[-idx]



# Granges must contain gene id in metadata column to create annotation files

genes@elementMetadata$id <- genes@elementMetadata@listData[["symbol"]]
genes_200up@elementMetadata$id <- genes_200up@elementMetadata@listData[["symbol"]]
genes_2000up@elementMetadata$id <- genes_2000up@elementMetadata@listData[["symbol"]]
genes_200up_gene@elementMetadata$id <- genes_200up_gene@elementMetadata@listData[["symbol"]]
genes_2000up_gene@elementMetadata$id <- genes_2000up_gene@elementMetadata@listData[["symbol"]]


## Create annotation file from Granges object for counting reads

annot_gene <- createAnnotationFile(genes)
annot_200up <- createAnnotationFile(genes_200up)
annot_2000up <- createAnnotationFile(genes_2000up)
annot_gene_200up <- createAnnotationFile(genes_200up_gene)
annot_gene_2000up <- createAnnotationFile(genes_2000up_gene)


# Make sure all annotation files are the same size

# use annot 2000 up as the base index

idx <- c(annot_2000up[["GeneID"]])

idx2 <- which(annot_200up[["GeneID"]] %in% idx)

annot_200up <- annot_200up[idx2,]
annot_gene <- annot_gene[idx2,]
annot_gene_200up <- annot_gene_200up[idx2,]



# run featureCounts to count up Tn5 insertions in each interval, for all replicates

snu_bam1 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S1_S6.lanemerge.dp.bam"
snu_bam2 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S2_S4.lanemerge.dp.bam"
snu_bam3 <- "/home/ubuntu1/atac/snu61/wt01/dp_bam/SNU61-WT-01-S3_S1.lanemerge.dp.bam"

## Gene ##
snu_r1_counts_gene <- featureCounts(
  files=snu_bam1,
  nthreads=20,
  annot.ext=annot_gene)

snu_r2_counts_gene <- featureCounts(
  files=snu_bam2,
  nthreads=20,
  annot.ext=annot_gene)

snu_r3_counts_gene <- featureCounts(
  files=snu_bam3,
  nthreads=20,
  annot.ext=annot_gene)

## 200 bp upstream
snu_r1_counts_200up <- featureCounts(
  files=snu_bam1,
  isPairedEnd=TRUE,
  nthreads=20,
  annot.ext=annot_200up)

snu_r2_counts_200up <- featureCounts(
  files=snu_bam2,
  nthreads=20,
  annot.ext=annot_200up)

snu_r3_counts_200up <- featureCounts(
  files=snu_bam3,
  nthreads=20,
  annot.ext=annot_200up)

## 2000 bp upstream
snu_r1_counts_2000up <- featureCounts(
  files=snu_bam1,
  isPairedEnd=TRUE,
  nthreads=20,
  annot.ext=annot_2000up)

snu_r2_counts_2000up <- featureCounts(
  files=snu_bam2,
  nthreads=20,
  annot.ext=annot_2000up)

snu_r3_counts_2000up <- featureCounts(
  files=snu_bam3,
  nthreads=20,
  annot.ext=annot_2000up)

## 200 up plus gene
snu_r1_counts_gene_200up <- featureCounts(
  files=snu_bam1,
  nthreads=20,
  annot.ext=annot_gene_200up)

snu_r2_counts_gene_200up <- featureCounts(
  files=snu_bam2,
  nthreads=20,
  annot.ext=annot_gene_200up)

snu_r3_counts_gene_200up <- featureCounts(
  files=snu_bam3,
  nthreads=20,
  annot.ext=annot_gene_200up)

## 2000 up and gene
snu_r1_counts_gene_2000up <- featureCounts(
  files=snu_bam1,
  nthreads=20,
  annot.ext=annot_gene_2000up)

snu_r2_counts_gene_2000up <- featureCounts(
  files=snu_bam2,
  nthreads=20,
  annot.ext=annot_gene_2000up)

snu_r3_counts_gene_2000up <- featureCounts(
  files=snu_bam3,
  nthreads=20,
  annot.ext=annot_gene_2000up)


##### LS1034 WT 01
ls1034_bam1 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S2_S3.lanemerge.dp.bam"
ls1034_bam2 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S1_S5.lanemerge.dp.bam"
ls1034_bam3 <- "/home/ubuntu1/atac/ls1034/wt01/dp_bam/LS1034-WT-01-S3_S2.lanemerge.dp.bam"

## Gene ##
ls1034_r1_counts_gene <- featureCounts(
  files=ls1034_bam1,
  nthreads=20,
  annot.ext=annot_gene)

ls1034_r2_counts_gene <- featureCounts(
  files=ls1034_bam2,
  nthreads=20,
  annot.ext=annot_gene)

ls1034_r3_counts_gene <- featureCounts(
  files=ls1034_bam3,
  nthreads=20,
  annot.ext=annot_gene)

## 200 bp upstream
ls1034_r1_counts_200up <- featureCounts(
  files=ls1034_bam1,
  nthreads=20,
  annot.ext=annot_200up)

ls1034_r2_counts_200up <- featureCounts(
  files=ls1034_bam2,
  nthreads=20,
  annot.ext=annot_200up)

ls1034_r3_counts_200up <- featureCounts(
  files=ls1034_bam3,
  nthreads=20,
  annot.ext=annot_200up)

## 2000 bp upstream
ls1034_r1_counts_2000up <- featureCounts(
  files=ls1034_bam1,
  nthreads=20,
  annot.ext=annot_2000up)

ls1034_r2_counts_2000up <- featureCounts(
  files=ls1034_bam2,
  nthreads=20,
  annot.ext=annot_2000up)

ls1034_r3_counts_2000up <- featureCounts(
  files=ls1034_bam3,
  nthreads=20,
  annot.ext=annot_2000up)

## 200 up plus gene
ls1034_r1_counts_gene_200up <- featureCounts(
  files=ls1034_bam1,
  nthreads=20,
  annot.ext=annot_gene_200up)

ls1034_r2_counts_gene_200up <- featureCounts(
  files=ls1034_bam2,
  nthreads=20,
  annot.ext=annot_gene_200up)

ls1034_r3_counts_gene_200up <- featureCounts(
  files=ls1034_bam3,
  nthreads=20,
  annot.ext=annot_gene_200up)

## 2000 up and gene
ls1034_r1_counts_gene_2000up <- featureCounts(
  files=ls1034_bam1,
  nthreads=20,
  annot.ext=annot_gene_2000up)

ls1034_r2_counts_gene_2000up <- featureCounts(
  files=ls1034_bam2,
  nthreads=20,
  annot.ext=annot_gene_2000up)

ls1034_r3_counts_gene_2000up <- featureCounts(
  files=ls1034_bam3,
  nthreads=20,
  annot.ext=annot_gene_2000up)

#### H508 WT 01
h508_bam1 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-2_S2.lanemerge.dp.bam"
h508_bam2 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-3_S1.lanemerge.dp.bam"
h508_bam3 <- "/home/ubuntu1/atac/h508/wt01/dp_bam/H508-1_S3.lanemerge.dp.bam"

## Gene ##
h508_r1_counts_gene <- featureCounts(
  files=h508_bam1,
  nthreads=20,
  annot.ext=annot_gene)

h508_r2_counts_gene <- featureCounts(
  files=h508_bam2,
  nthreads=20,
  annot.ext=annot_gene)

h508_r3_counts_gene <- featureCounts(
  files=h508_bam3,
  nthreads=20,
  annot.ext=annot_gene)

## 200 bp upstream
h508_r1_counts_200up <- featureCounts(
  files=h508_bam1,
  nthreads=20,
  annot.ext=annot_200up)

h508_r2_counts_200up <- featureCounts(
  files=h508_bam2,
  nthreads=20,
  annot.ext=annot_200up)

h508_r3_counts_200up <- featureCounts(
  files=h508_bam3,
  nthreads=20,
  annot.ext=annot_200up)

## 2000 bp upstream
h508_r1_counts_2000up <- featureCounts(
  files=h508_bam1,
  nthreads=20,
  annot.ext=annot_2000up)

h508_r2_counts_2000up <- featureCounts(
  files=h508_bam2,
  nthreads=20,
  annot.ext=annot_2000up)

h508_r3_counts_2000up <- featureCounts(
  files=h508_bam3,
  nthreads=20,
  annot.ext=annot_2000up)

## 200 up plus gene
h508_r1_counts_gene_200up <- featureCounts(
  files=h508_bam1,
  nthreads=20,
  annot.ext=annot_gene_200up)

h508_r2_counts_gene_200up <- featureCounts(
  files=h508_bam2,
  nthreads=20,
  annot.ext=annot_gene_200up)

h508_r3_counts_gene_200up <- featureCounts(
  files=h508_bam3,
  nthreads=20,
  annot.ext=annot_gene_200up)

## 2000 up and gene
h508_r1_counts_gene_2000up <- featureCounts(
  files=h508_bam1,
  nthreads=20,
  annot.ext=annot_gene_2000up)

h508_r2_counts_gene_2000up <- featureCounts(
  files=h508_bam2,
  nthreads=20,
  annot.ext=annot_gene_2000up)

h508_r3_counts_gene_2000up <- featureCounts(
  files=h508_bam3,
  nthreads=20,
  annot.ext=annot_gene_2000up)


## Copy the data to a matrix


one <- paste0("snu61",c("gene_r1","gene_r2","gene_r3","200up_r1","200up_r2","200up_r3","2000up_r1","2000up_r2","2000up_r3","gene+200up_r1","gene+200up_r2","gene+200up_r3","gene+2000up_r1","gene+2000up_r2","gene+2000up_r3"))
two <- paste0("ls1034",c("gene_r1","gene_r2","gene_r3","200up_r1","200up_r2","200up_r3","2000up_r1","2000up_r2","2000up_r3","gene+200up_r1","gene+200up_r2","gene+200up_r3","gene+2000up_r1","gene+2000up_r2","gene+2000up_r3"))
three <- paste0("h508",c("gene_r1","gene_r2","gene_r3","200up_r1","200up_r2","200up_r3","2000up_r1","2000up_r2","2000up_r3","gene+200up_r1","gene+200up_r2","gene+200up_r3","gene+2000up_r1","gene+2000up_r2","gene+2000up_r3"))

four <- paste0(c("snu61_", "ls1034_", "h508_"), "gene_avg")
five <- paste0(c("snu61_", "ls1034_", "h508_"), "200bp_avg")
six <- paste0(c("snu61_", "ls1034_", "h508_"), "2000bp_avg")
seven <- paste0(c("snu61_", "ls1034_", "h508_"), "gene+200_avg")
eight <- paste0(c("snu61_", "ls1034_", "h508_"), "gene+2000_avg")

nine <- c("gene_diff_h508_ls1034","gene_diff_h508_snu61","gene_diff_snu61_ls1034")
ten <- c("200_diff_h508_ls1034","200_diff_h508_snu61","200_diff_snu61_ls1034")
eleven <- c("2000_diff_h508_ls1034","2000_diff_h508_snu61","2000_diff_snu61_ls1034")
twelve <- c("gene+200_diff_h508_ls1034","gene+200_diff_h508_snu61","gene+200_diff_snu61_ls1034")
thirteen <- c("gene+2000_diff_h508_ls1034","gene+2000_diff_h508_snu61","gene+2000_diff_snu61_ls1034")

colnames <- c(one,two,three,four,five,six,seven,eight,nine,ten,eleven,twelve,thirteen)
rm(one,two,three,four,five,six,seven,eight,nine,ten,eleven,twelve,thirteen)

# make the matrix
items <- length(annot_gene[,1])
#
signal_matrix <- matrix(data = NA, items, 75)
#
colnames(signal_matrix) <- colnames
# apply rownames to matrix
rownames(signal_matrix) <- c(annot_gene[["GeneID"]])


# Copy signals to matrix

for (c in 1:items){
  
  ## SNU61 raw data
  signal_matrix[[c,1]] <- as.integer(snu_r1_counts_gene[["counts"]][[c]])
  signal_matrix[[c,2]] <- as.integer(snu_r2_counts_gene[["counts"]][[c]])
  signal_matrix[[c,3]] <- as.integer(snu_r3_counts_gene[["counts"]][[c]])
  signal_matrix[[c,4]] <- as.integer(snu_r1_counts_200up[["counts"]][[c]])
  signal_matrix[[c,5]] <- as.integer(snu_r2_counts_200up[["counts"]][[c]])
  signal_matrix[[c,6]] <- as.integer(snu_r3_counts_200up[["counts"]][[c]])
  signal_matrix[[c,7]] <- as.integer(snu_r1_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,8]] <- as.integer(snu_r2_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,9]] <- as.integer(snu_r3_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,10]] <- as.integer(snu_r1_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,11]] <- as.integer(snu_r2_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,12]] <- as.integer(snu_r3_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,13]] <- as.integer(snu_r1_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,14]] <- as.integer(snu_r2_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,15]] <- as.integer(snu_r3_counts_gene_2000up[["counts"]][[c]])
  
  ## LS1034 raw data
  signal_matrix[[c,16]] <- as.integer(ls1034_r1_counts_gene[["counts"]][[c]])
  signal_matrix[[c,17]] <- as.integer(ls1034_r2_counts_gene[["counts"]][[c]])
  signal_matrix[[c,18]] <- as.integer(ls1034_r3_counts_gene[["counts"]][[c]])
  signal_matrix[[c,19]] <- as.integer(ls1034_r1_counts_200up[["counts"]][[c]])
  signal_matrix[[c,20]] <- as.integer(ls1034_r2_counts_200up[["counts"]][[c]])
  signal_matrix[[c,21]] <- as.integer(ls1034_r3_counts_200up[["counts"]][[c]])
  signal_matrix[[c,22]] <- as.integer(ls1034_r1_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,23]] <- as.integer(ls1034_r2_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,24]] <- as.integer(ls1034_r3_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,25]] <- as.integer(ls1034_r1_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,26]] <- as.integer(ls1034_r2_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,27]] <- as.integer(ls1034_r3_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,28]] <- as.integer(ls1034_r1_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,29]] <- as.integer(ls1034_r2_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,30]] <- as.integer(ls1034_r3_counts_gene_2000up[["counts"]][[c]])
  
  ## H508 raw data
  signal_matrix[[c,31]] <- as.integer(h508_r1_counts_gene[["counts"]][[c]])
  signal_matrix[[c,32]] <- as.integer(h508_r2_counts_gene[["counts"]][[c]])
  signal_matrix[[c,33]] <- as.integer(h508_r3_counts_gene[["counts"]][[c]])
  signal_matrix[[c,34]] <- as.integer(h508_r1_counts_200up[["counts"]][[c]])
  signal_matrix[[c,35]] <- as.integer(h508_r2_counts_200up[["counts"]][[c]])
  signal_matrix[[c,36]] <- as.integer(h508_r3_counts_200up[["counts"]][[c]])
  signal_matrix[[c,37]] <- as.integer(h508_r1_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,38]] <- as.integer(h508_r2_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,39]] <- as.integer(h508_r3_counts_2000up[["counts"]][[c]])
  signal_matrix[[c,40]] <- as.integer(h508_r1_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,41]] <- as.integer(h508_r2_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,42]] <- as.integer(h508_r3_counts_gene_200up[["counts"]][[c]])
  signal_matrix[[c,43]] <- as.integer(h508_r1_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,44]] <- as.integer(h508_r2_counts_gene_2000up[["counts"]][[c]])
  signal_matrix[[c,45]] <- as.integer(h508_r3_counts_gene_2000up[["counts"]][[c]])
  
  ## Gene signal averages (snu61, ls1034, h508)
  signal_matrix[[c,46]] <- as.integer(mean(signal_matrix[c,1],signal_matrix[c,2],signal_matrix[c,3]))
  signal_matrix[[c,47]] <- as.integer(mean(signal_matrix[c,16],signal_matrix[c,17],signal_matrix[c,18]))
  signal_matrix[[c,48]] <- as.integer(mean(signal_matrix[c,31],signal_matrix[c,32],signal_matrix[c,33]))
  
  ## 200 bp up signal averages (snu61, ls1034, h508)
  signal_matrix[[c,49]] <- as.integer(mean(signal_matrix[c,4],signal_matrix[c,5],signal_matrix[c,6]))
  signal_matrix[[c,50]] <- as.integer(mean(signal_matrix[c,19],signal_matrix[c,20],signal_matrix[c,21]))
  signal_matrix[[c,51]] <- as.integer(mean(signal_matrix[c,34],signal_matrix[c,35],signal_matrix[c,36]))
  
  ## 2000 bp up signal averages (snu61, ls1034, h508)
  signal_matrix[[c,52]] <- as.integer(mean(signal_matrix[c,7],signal_matrix[c,8],signal_matrix[c,9]))
  signal_matrix[[c,53]] <- as.integer(mean(signal_matrix[c,22],signal_matrix[c,23],signal_matrix[c,24]))
  signal_matrix[[c,54]] <- as.integer(mean(signal_matrix[c,37],signal_matrix[c,38],signal_matrix[c,39]))
  
  ## 200 bp up and gene signal averages (snu61, ls1034, h508)
  signal_matrix[[c,55]] <- as.integer(mean(signal_matrix[c,10],signal_matrix[c,11],signal_matrix[c,12]))
  signal_matrix[[c,56]] <- as.integer(mean(signal_matrix[c,25],signal_matrix[c,26],signal_matrix[c,27]))
  signal_matrix[[c,57]] <- as.integer(mean(signal_matrix[c,40],signal_matrix[c,41],signal_matrix[c,42]))
  
  ## 2000 bp and up gene signal averages (snu61, ls1034, h508)
  signal_matrix[[c,58]] <- as.integer(mean(signal_matrix[c,13],signal_matrix[c,14],signal_matrix[c,15]))
  signal_matrix[[c,59]] <- as.integer(mean(signal_matrix[c,28],signal_matrix[c,29],signal_matrix[c,30]))
  signal_matrix[[c,60]] <- as.integer(mean(signal_matrix[c,43],signal_matrix[c,44],signal_matrix[c,45]))
  
  #### Pairwise differences
  
  ## gene (h508/ls1034, h508/snu61, snu61/ls1034)
  signal_matrix[[c,61]] <- abs(signal_matrix[[c,48]]-signal_matrix[[c,47]])
  signal_matrix[[c,62]] <- abs(signal_matrix[[c,48]]-signal_matrix[[c,46]])
  signal_matrix[[c,63]] <- abs(signal_matrix[[c,47]]-signal_matrix[[c,46]])
  
  ## 200 bp upstream (h508/ls1034, h508/snu61, snu61/ls1034)
  signal_matrix[[c,64]] <- abs(signal_matrix[[c,51]]-signal_matrix[[c,50]])
  signal_matrix[[c,65]] <- abs(signal_matrix[[c,51]]-signal_matrix[[c,49]])
  signal_matrix[[c,66]] <- abs(signal_matrix[[c,50]]-signal_matrix[[c,49]])
  
  ## 2000 bp upstream (h508/ls1034, h508/snu61, snu61/ls1034)
  signal_matrix[[c,67]] <- abs(signal_matrix[[c,54]]-signal_matrix[[c,53]])
  signal_matrix[[c,68]] <- abs(signal_matrix[[c,54]]-signal_matrix[[c,52]])
  signal_matrix[[c,69]] <- abs(signal_matrix[[c,53]]-signal_matrix[[c,52]])
  
  ## 200 bp upstream and gene (h508/ls1034, h508/snu61, snu61/ls1034)
  signal_matrix[[c,70]] <- abs(signal_matrix[[c,57]]-signal_matrix[[c,56]])
  signal_matrix[[c,71]] <- abs(signal_matrix[[c,57]]-signal_matrix[[c,55]])
  signal_matrix[[c,72]] <- abs(signal_matrix[[c,56]]-signal_matrix[[c,55]])
  
  ## 2000 bp upstream and gene (h508/ls1034, h508/snu61, snu61/ls1034)
  signal_matrix[[c,73]] <- abs(signal_matrix[[c,60]]-signal_matrix[[c,59]])
  signal_matrix[[c,74]] <- abs(signal_matrix[[c,60]]-signal_matrix[[c,58]])
  signal_matrix[[c,75]] <- abs(signal_matrix[[c,59]]-signal_matrix[[c,58]])
}


# Convert values of interest to dataframe to make descriptive plots

# snu61, ls1034, h508
df <- data.frame(
  "snu61_gene"=signal_matrix[,46],
  "ls1034_gene"=signal_matrix[,47],
  "h508_gene"=signal_matrix[,48],
  #
  "snu61_200"=signal_matrix[,49],
  "ls1034_200"=signal_matrix[,50],
  "h508_200"=signal_matrix[,51],
  #
  "snu61_2000"=signal_matrix[,52],
  "ls1034_2000"=signal_matrix[,53],
  "h508_2000"=signal_matrix[,54],
  #
  "snu61_200_gene"=signal_matrix[,55],
  "ls1034_200_gene"=signal_matrix[,56],
  "h508_200_gene"=signal_matrix[,57],
  #
  "snu61_2000_gene"=signal_matrix[,58],
  "ls1034_2000_gene"=signal_matrix[,59],
  "h508_2000_gene"=signal_matrix[,60]
)



# Make descriptive statistics and plots from data

options(scipen=999)  # turn-off scientific notation like 1e+4


## calculate quantiles ###
quantiles <- list()
# SNU61 quantiles
quantiles$snu61_gene <- quantile(df[["snu61_gene"]], na.rm=TRUE)
quantiles$snu61_200 <- quantile(df[["snu61_200"]], na.rm=TRUE)
quantiles$snu61_2000 <- quantile(df[["snu61_2000"]], na.rm=TRUE)
quantiles$snu61_200_gene <- quantile(df[["snu61_200_gene"]], na.rm=TRUE)
quantiles$snu61_2000_gene <- quantile(df[["snu61_2000_gene"]], na.rm=TRUE)
# LS1034 quantiles
quantiles$ls1034_gene <- quantile(df[["ls1034_gene"]], na.rm=TRUE)
quantiles$ls1034_200 <- quantile(df[["ls1034_200"]], na.rm=TRUE)
quantiles$ls1034_2000 <- quantile(df[["ls1034_2000"]], na.rm=TRUE)
quantiles$ls1034_200_gene <- quantile(df[["ls1034_200_gene"]], na.rm=TRUE)
quantiles$ls1034_2000_gene <- quantile(df[["ls1034_2000_gene"]], na.rm=TRUE)
# H508 quantiles
quantiles$h508_gene <- quantile(df[["h508_gene"]], na.rm=TRUE)
quantiles$h508_200 <- quantile(df[["h508_200"]], na.rm=TRUE)
quantiles$h508_2000 <- quantile(df[["h508_2000"]], na.rm=TRUE)
quantiles$h508_200_gene <- quantile(df[["h508_200_gene"]], na.rm=TRUE)
quantiles$h508_2000_gene <- quantile(df[["h508_2000_gene"]], na.rm=TRUE)


### make basic scatter plots to show range of signal values ###
# SNU61
plot.default(x=df[["snu61_gene"]], y=NULL, ylab="READS GENE", xlab="SNU-61")
plot.default(x=df[["snu61_200"]], y=NULL, ylab="READS 200", xlab="SNU-61")
plot.default(x=df[["snu61_2000"]], y=NULL, ylab="READS 2000", xlab="SNU-61")
plot.default(x=df[["snu61_200_gene"]], y=NULL, ylab="READS 200 GENE", xlab="SNU-61")
plot.default(x=df[["snu61_2000_gene"]], y=NULL, ylab="READS 2000 GENE", xlab="SNU-61")
# LS1034
plot.default(x=df[["ls1034_gene"]], y=NULL, ylab="READS GENE", xlab="LS1034")
plot.default(x=df[["ls1034_200"]], y=NULL, ylab="READS 200", xlab="LS1034")
plot.default(x=df[["ls1034_2000"]], y=NULL, ylab="READS 2000", xlab="LS1034")
plot.default(x=df[["ls1034_200_gene"]], y=NULL, ylab="READS 200 GENE", xlab="LS1034")
plot.default(x=df[["ls1034_2000_gene"]], y=NULL, ylab="READS 2000 GENE", xlab="LS1034")
# H508
plot.default(x=df[["h508_gene"]], y=NULL, ylab="READS GENE", xlab="H508")
plot.default(x=df[["h508_200"]], y=NULL, ylab="READS 200", xlab="H508")
plot.default(x=df[["h508_2000"]], y=NULL, ylab="READS 2000", xlab="H508")
plot.default(x=df[["h508_200_gene"]], y=NULL, ylab="READS 200 GENE", xlab="H508")
plot.default(x=df[["h508_2000_gene"]], y=NULL, ylab="READS 2000 GENE", xlab="H508")


### Make scatterplots comparing two columns (categories or cell lines) ###
# Scatterplot
theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(df, aes(snu61_gene, snu61_200))
g + geom_jitter(width = .5, size=1) +
  labs(y="200 bp upstream", 
       x="gene", 
       title="Read counts gene v 200 bp upstream")

#** log transforms? **#

### Make histograms ###
hist(df[["snu61_gene"]], breaks = 50, col="blue", main="snu61 gene")
hist(df[["ls1034_gene"]], breaks = 50, col="blue", main="ls1034 gene")
hist(df[["h508_gene"]], breaks = 50, col="blue", main="h508 gene")

hist(df[["snu61_200"]], breaks = 50, col="blue", main="snu61 200")

hist(df[["snu61_2000"]], breaks = 50, col="blue", main="snu61 2000")

hist(df[["snu61_200_gene"]], breaks = 50, col="blue", main="snu61 200 gene")

hist(df[["snu61_2000_gene"]], breaks = 50, col="blue", main="snu61 2000 gene")




## archived code (for now) subset granges for specific targets only (like from a text file) **

textfile <- "/home/rstudio1/atac/h508/smad4_entrez_targets.txt"
con <- file(description=textfile, open="r")

## Copy list of ENTREZ target IDS from text file
smad4_targets <- c()
for (i in 1:7897){
  smad4_targets[i] <- scan(file=con, nlines=1, quiet=TRUE)
}

#### end ####
close.connection(con=con)
#### Subset the granges object for specific targets ####
targets <- c(which(g@elementMetadata@listData[["entrezid"]] %in% smad4_targets))
g_smad4 <- g
# this metadata column is required to make the annot object
g_smad4@elementMetadata$id <- g_smad4@elementMetadata@listData[["gene_name"]]
#### end ####


reg_200up_counts <- list()
#
reg_200up_counts$snu61_r1 <- snu_r1_counts_200up
reg_200up_counts$snu61_r2 <- snu_r2_counts_200up
reg_200up_counts$snu61_r3 <- snu_r3_counts_200up
reg_200up_counts$ls1034_r1 <- ls1034_r1_counts_200up
reg_200up_counts$ls1034_r2 <- ls1034_r2_counts_200up
reg_200up_counts$ls1034_r3 <- ls1034_r3_counts_200up
reg_200up_counts$h508_r1 <- h508_r1_counts_200up
reg_200up_counts$h508_r2 <- h508_r2_counts_200up
reg_200up_counts$h508_r3 <- h508_r3_counts_200up


reg_2000up_counts <- list()
#
reg_2000up_counts$snu61_r1 <- snu_r1_counts_2000up
reg_2000up_counts$snu61_r2 <- snu_r2_counts_2000up
reg_2000up_counts$snu61_r3 <- snu_r3_counts_2000up
reg_2000up_counts$ls1034_r1 <- ls1034_r1_counts_2000up
reg_2000up_counts$ls1034_r2 <- ls1034_r2_counts_2000up
reg_2000up_counts$ls1034_r3 <- ls1034_r3_counts_2000up
reg_2000up_counts$h508_r1 <- h508_r1_counts_2000up
reg_2000up_counts$h508_r2 <- h508_r2_counts_2000up
reg_2000up_counts$h508_r3 <- h508_r3_counts_2000up


num_genes <- length(reg_200up_counts[["snu61_r1"]][["counts"]])
gene_names <- reg_200up_counts[["snu61_r1"]][["annotation"]][["GeneID"]]
#
reg_200_norm_counts <- list()
reg_200_norm_counts$snu61 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(reg_200_norm_counts$snu61) <- c("rep1","rep2","rep3","avg","total")
rownames(reg_200_norm_counts$snu61) <- c(gene_names)
reg_200_norm_counts$ls1034 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(reg_200_norm_counts$ls1034) <- c("rep1","rep2","rep3","avg","total")
rownames(reg_200_norm_counts$ls1034) <- c(gene_names)
reg_200_norm_counts$h508 <- matrix(data=NA, ncol=5, nrow=num_genes)
colnames(reg_200_norm_counts$h508) <- c("rep1","rep2","rep3","avg","total")
rownames(reg_200_norm_counts$h508) <- c(gene_names)
#
for (a in 1:num_genes){
  #
  reg_200_norm_counts$snu61[a,1] <- as.integer(reg_200up_counts[["snu61_r1"]][["counts"]][[a]]/replicate_normalization[1,2])
  reg_200_norm_counts$snu61[a,2] <- as.integer(reg_200up_counts[["snu61_r2"]][["counts"]][[a]]/replicate_normalization[2,2])
  reg_200_norm_counts$snu61[a,3] <- as.integer(reg_200up_counts[["snu61_r3"]][["counts"]][[a]]/replicate_normalization[3,2])
  reg_200_norm_counts$snu61[a,4] <- ((reg_200_norm_counts$snu61[a,1]+reg_200_norm_counts$snu61[a,2]+reg_200_norm_counts$snu61[a,3])/3)
  reg_200_norm_counts$snu61[a,5] <- reg_200_norm_counts$snu61[a,1]+reg_200_norm_counts$snu61[a,2]+reg_200_norm_counts$snu61[a,3]
  #
  reg_200_norm_counts$ls1034[a,1] <- as.integer(reg_200up_counts[["ls1034_r1"]][["counts"]][[a]]/replicate_normalization[4,2])
  reg_200_norm_counts$ls1034[a,2] <- as.integer(reg_200up_counts[["ls1034_r2"]][["counts"]][[a]]/replicate_normalization[5,2])
  reg_200_norm_counts$ls1034[a,3] <- as.integer(reg_200up_counts[["ls1034_r3"]][["counts"]][[a]]/replicate_normalization[6,2])
  reg_200_norm_counts$ls1034[a,4] <- ((reg_200_norm_counts$ls1034[a,1]+reg_200_norm_counts$ls1034[a,2]+reg_200_norm_counts$ls1034[a,3])/3)
  reg_200_norm_counts$ls1034[a,5] <- reg_200_norm_counts$ls1034[a,1]+reg_200_norm_counts$ls1034[a,2]+reg_200_norm_counts$ls1034[a,3]
  #
  reg_200_norm_counts$h508[a,1] <- as.integer(reg_200up_counts[["h508_r1"]][["counts"]][[a]]/replicate_normalization[7,2])
  reg_200_norm_counts$h508[a,2] <- as.integer(reg_200up_counts[["h508_r2"]][["counts"]][[a]]/replicate_normalization[8,2])
  reg_200_norm_counts$h508[a,3] <- as.integer(reg_200up_counts[["h508_r3"]][["counts"]][[a]]/replicate_normalization[9,2])
  reg_200_norm_counts$h508[a,4] <- ((reg_200_norm_counts$h508[a,1]+reg_200_norm_counts$h508[a,2]+reg_200_norm_counts$h508[a,3])/3)
  reg_200_norm_counts$h508[a,5] <- reg_200_norm_counts$h508[a,1]+reg_200_norm_counts$h508[a,2]+reg_200_norm_counts$h508[a,3]}
#

#
reg200_signal_matrix <- matrix(data=NA,nrow=num_genes,ncol=5)
rownames(reg200_signal_matrix) <- c(gene_names)
colnames(reg200_signal_matrix) <- c("snu61","ls1034","h508","var","sd")
#
for (b in 1:num_genes){
  #
  reg200_signal_matrix[b,1] <- as.integer((reg_200_norm_counts[["snu61"]][b,4]/sample_normalization[1,2]))
  reg200_signal_matrix[b,2] <- as.integer((reg_200_norm_counts[["ls1034"]][b,4]/sample_normalization[2,2]))
  reg200_signal_matrix[b,3] <- as.integer((reg_200_norm_counts[["h508"]][b,4]/sample_normalization[3,2]))
  reg200_signal_matrix[b,4] <- stats::var(reg200_signal_matrix[b,1:3])
  reg200_signal_matrix[b,5] <- stats::sd(reg200_signal_matrix[b,1:3])
}
#
coad_mr <- c("TCF7","MNX1","POU5F1B","ESRRA","CDX2","HNF4A","GMEB2","HOXA3","OVOL1","ASCL2","ZSWIM1","CBFA2T2","TARBP1","KLHL31","GTF2IRD1","ZNF696","NUFIP1","SUPT20H","KAT8","KAT2A","TAF4","RBM39","ZMYND8","L3MBTL1","PLAGL2","NCOA5","ZSWIM3","ADNP","ASXL1")
ind <- which(rownames(reg200_signal_matrix) %in% coad_mr)
reg200_mr_counts <- reg200_signal_matrix[ind,]


# Violin plots
# prep data
violin_data_reg200 <- data.frame("Signal" = c(reg200_mr_counts[,1],reg200_mr_counts[,2],reg200_mr_counts[,3]),
                                 "Sample" = c(rep("snu61", times=29), rep("ls1034",times=29),rep("h508",times=29)))
#
violin_plot_reg200 <- ggplot(violin_data_reg200, aes(x=Sample, y=Signal, fill=Sample)) + 
  geom_violin()
violin_plot_reg200


heatmap_data_reg200 <- data.frame("Signal" = c(reg200_mr_counts[,1],reg200_mr_counts[,2],reg200_mr_counts[,3]),
                                  "Sample" = c(rep("snu61", times=29), rep("ls1034",times=29),rep("h508",times=29)),
                                  "Gene" = rep(mr_names))

heatmap.plot.reg200 <- ggplot(data = heatmap_data_reg200, aes(x = Sample, y = Gene)) +
  geom_tile(aes(fill = Signal)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  theme(axis.text.y = element_text(size = 6))

# Preview the heatmap
print(heatmap.plot.reg200)

#
mr_all_reps_reg200 <- matrix(data=NA, nrow=29, ncol=9)
rownames(mr_all_reps_reg200) <- mr_names
colnames(mr_all_reps_reg200) <- c("snu61_r1","snu61_r2","snu61_r3","ls1034_r1","ls1034_r2","ls1034_r3","h508_r1","h508_r2","h508_r3")
#
mr_all_reps_reg200[,1] <- reg_200_norm_counts[["snu61"]][ind,1]
mr_all_reps_reg200[,2] <- reg_200_norm_counts[["snu61"]][ind,2]
mr_all_reps_reg200[,3] <- reg_200_norm_counts[["snu61"]][ind,3]
mr_all_reps_reg200[,4] <- reg_200_norm_counts[["ls1034"]][ind,1]
mr_all_reps_reg200[,5] <- reg_200_norm_counts[["ls1034"]][ind,2]
mr_all_reps_reg200[,6] <- reg_200_norm_counts[["ls1034"]][ind,3]
mr_all_reps_reg200[,7] <- reg_200_norm_counts[["h508"]][ind,1]
mr_all_reps_reg200[,8] <- reg_200_norm_counts[["h508"]][ind,2]
mr_all_reps_reg200[,9] <- reg_200_norm_counts[["h508"]][ind,3]
#
heatmap_data_reps_reg200 <- data.frame(
  "Signal" = c(mr_all_reps_reg200[,1],mr_all_reps_reg200[,2],mr_all_reps_reg200[,3],mr_all_reps_reg200[,4],mr_all_reps_reg200[,5],mr_all_reps_reg200[,6],mr_all_reps_reg200[,7],mr_all_reps_reg200[,8],mr_all_reps_reg200[,9]),
  "Sample" = c(rep("snu61_r1", times=29), rep("snu61_r2", times=29), rep("snu61_r3", times=29), rep("ls1034_r1",times=29), rep("ls1034_r2",times=29), rep("ls1034_r3",times=29), rep("h508_r1",times=29), rep("h508_r2",times=29), rep("h508_r3",times=29)),
  "Gene" = rep(mr_names))
#
heatmap.plot.reps.reg200 <- ggplot(data = heatmap_data_reps_reg200, aes(x = Sample, y = Gene)) +
  geom_tile(aes(fill = Signal)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  theme(axis.text.y = element_text(size = 6))

# Preview the heatmap
print(heatmap.plot.reps.reg200)





## load libs
library(ChIPpeakAnno)
library(genomation)
library(org.Hs.eg.db)

## import
snu61_peak <- readBed("/home/ubuntu1/atac/snu61/wt01/preprocessing/13allpeaks/SNU61-WT-01.all_summits.bed")
snu61_peak2 <- readBed("/home/ubuntu1/atac/snu61/wt01/preprocessing/13allpeaks/SNU61-WT-01.all_peaks.narrowPeak")




##
snu61.anno <- annotatePeakInBatch(snu61_peak, AnnotationData=TSS.human.GRCh38)

## Add gene symbols
snu61.anno <- addGeneIDs(annotatedPeak=snu61.anno, 
                         orgAnn="org.Hs.eg.db", 
                         IDs2Add="symbol")

## Find unique
snu61.genes <- unique(snu61.anno@elementMetadata@listData[["symbol"]])


x <- keepStandardChromosomes(snu61_peak)


## Write to file
write.table(snu61_ls1034_genes, file = "/home/ubuntu1/atac/peaks/snu61_ls1034_genes.txt", quote = FALSE, sep = ",", eol = "\n", row.names = FALSE, col.names = FALSE)

## Overlaps
ov <- findOverlaps(snu61_peak2, TSS.human.GRCh38)

######################################################################################

## annotate the peaks with precompiled ensembl annotation
data(TSS.human.GRCh38)
## Do some coercion to convert the TSS.human GRANGEs seqnames to standard format (with chr)
newseqs <- paste0("chr", c(1:22, "X", "Y"))
names(newseqs) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                    "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                    "21", "22", "X", "Y")
## Rename the sequences
newTSS <- renameSeqlevels(TSS.human.GRCh38, newseqs)
## Find overlaps
nov <- findOverlaps(snu61_peak, newTSS)
nov2 <- unique(nov@to)
##
genes <- TSS.human.GRCh38@elementMetadata@listData[["description"]][nov2]
genes <- genes[-which(is.null(genes) == FALSE)]
genes <- genes[-which(genes == "")]

