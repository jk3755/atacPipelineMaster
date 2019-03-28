## Install/load required libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("stringi", suppressUpdates = TRUE)
#biocLite("GenomicRanges", suppressUpdates = TRUE)
#biocLite("ensembldb", suppressUpdates = TRUE)
#biocLite("EnsDb.Hsapiens.v86", suppressUpdates = TRUE)
#biocLite("S4Vectors", suppressUpdates = TRUE)
#biocLite("Rsubread", suppressUpdates = TRUE)
#biocLite("Repitools", suppressUpdates = TRUE)
#biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene", suppressUpdates = TRUE)
#biocLite("genomation", suppressUpdates = TRUE)
#install.packages("ggplot2")
##
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(genomation)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(S4Vectors)
library(Rsubread)
library(Repitools)
library(ggplot2)
library(crayon)
library(Rsamtools)


## Shorten variable name for TxDb database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#### The input files
## First file is provided by Jeremy, contains RNAseq data
inputFile <- "/home/ubuntu1/atac/jeremy/snu61_tfc"
expData <- read.table(inputFile, header=TRUE, sep = ',')
## Second file is called peaks from ATACseq data, in .bed format
bedFile <- "/home/ubuntu1/atac/snu61/wt01/preprocessing/13allpeaks/SNU61-WT-01.all_peaks.narrowPeak" # is it better to use narrowPeak or summits here?
snu61Peaks <- readBed(bedFile, track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
## Keep standard ranges only
snu61Peaks <- keepStandardChromosomes(snu61Peaks, pruning.mode="coarse")
## The input bam files, for calculating raw signals
inputBam <- "/home/ubuntu1/atac/snu61/wt01/preprocessing/12all/SNU61-WT-01.all.bam"


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

