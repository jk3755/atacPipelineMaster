##
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#BiocManager::install("viper")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("mygene")
#install.packages("rlist")

##
cat("Loading libraries...", "\n")
suppressMessages(library(aracne.networks))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(annotate))
suppressMessages(library(rlist))
suppressMessages(library(viper))
suppressMessages(library(GenomicRanges))
suppressMessages(library(mygene))

##
cat("Setting snakemake vars...", "\n")
inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
entrezid <- snakemake@wildcards[["entrez"]]

## load file
load(inputfile)


# ESRRA - 2101
################### ENTER ENTREZ ID AND RUN TO GET TARGETS ##############################
## get the ARACNe interactome
coad_interactome <- aracne.networks::reguloncoad
## get the targets
com <- paste0("targets <- names(coad_interactome[['", entrezid, "']][['tfmode']])")
eval(parse(text = com))
## Retrieve info for the network targets
target_list <- list()
for (a in 1:length(targets)){
  target_list[a] <- mygene::getGene(geneid = targets[a], fields = "all")}
## count the number of gene locations we have
loc_count <- 0
for (a in 1:length(targets)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      loc_count <- (loc_count + length(target_list[[a]][["genomic_pos"]]))
    } else {
      loc_count <- (loc_count + 1)}}}
## Retrieve the genomic coordinates of all network targets
target_locations <- matrix(data = NA, nrow = loc_count, ncol = 4)
colnames(target_locations) <- c("gene", "chr", "start", "end")
idx <- 1
#
for (a in 1:length(target_list)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      for (b in 1:length(target_list[[a]][["genomic_pos"]])){
        target_locations[idx,1] <- target_list[[a]][["symbol"]]
        target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][[b]][["chr"]])
        target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][[b]][["start"]]
        target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][[b]][["end"]]
        idx <- (idx+1)
      } # end for (b in 1:length(target_list[[a]][["genomic_pos"]]))
    } else {
      target_locations[idx,1] <- target_list[[a]][["symbol"]]
      target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][["chr"]])
      target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][["start"]]
      target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][["end"]]
      idx <- (idx+1)}}}
## Convert the genomic locations of the targets into GRanges
gr <- GRanges(
  seqnames = target_locations[,2],
  ranges = IRanges(start = as.numeric(target_locations[,3]), end = as.numeric(target_locations[,4])))
#strand = c(rep("+", times = num_names)),
#seqinfo = Seqinfo(genome="hg38"),
#score = c(rep(1, times = num_names)))
## prune to standard xsomes
gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
##
gr1 <- gr
gr2 <- gr
##
start(gr1) <- (start(gr1) - 2000)
width(gr1) <- (width(gr1) + 2000)
##
start(gr2) <- (start(gr2) - 50000)
width(gr2) <- (width(gr2) + 50000)
#################################################################################################################

#### PARSED SITES ####
parsed_overlap1 <- length(info[["overlap1"]]@from)
parsed_overlap2 <- length(info[["overlap1"]]@from)

## get binding sites of merged signals
sites <- mergedMotifs[["sites"]]


## Find overlaps
## YOU WILL WANT TO FIND THE OVERLAPS FOR A RANGE OF EXTENDED VALUES PAST THE TARGETS, GRAPH IT
overlaps1 <- findOverlaps(gr1, sites)
overlaps2 <- findOverlaps(gr2, sites)

info <- list()
info$overlap1 <- overlaps1
info$overlap2 <- overlaps2

save(info, file = outputfile)


## Fishers test
o1 <- length(info[["overlap1"]]@from)
o2 <- length(info[["overlap2"]]@from)
numtargets <- info[["overlap1"]]@nLnode
numsites <- info[["overlap1"]]@nRnode

percent1 <- (o1/numsites*100)

#
motif1sites <- parsedSitesInfo[["bfPassPeakSites"]]
overlap3 <- findOverlaps(gr1, sites)
o3 <- length(overlap3@from)
overlap4 <- findOverlaps(gr2, sites)
o4 <- length(overlap4@from)









