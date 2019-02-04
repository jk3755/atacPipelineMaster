##
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#BiocManager::install("viper")
#BiocManager::install("AnnotationDbi")
#install.packages("rlist")
##
cat("Loading libraries...", "\n")
library(aracne.networks)
library(org.Hs.eg.db)
library(annotate)
library(rlist)
library(viper)
library(GenomicRanges)


##
cat("Setting snakemake vars...", "\n")
inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]



## get the ARACNe interactome
coad_interactome <- aracne.networks::reguloncoad


## get the mnx1 targets
cdx2_targets <- names(coad_interactome[["1045"]][["tfmode"]])


## Retrieve info for the network targets
target_list <- list()
for (a in 1:length(cdx2_targets)){
  target_list[a] <- mygene::getGene(geneid = cdx2_targets[a], fields = "all")}


## count the number of gene locations we have
loc_count <- 0
for (a in 1:length(cdx2_targets)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      loc_count <- (loc_count + length(target_list[[a]][["genomic_pos"]]))
    } else {
      loc_count <- (loc_count + 1)}}}


## Retrieve the genomic coordinates of all network targets
## retrieve all genomic coords here, can later trim non standard ones more easily with GRanges functions
target_locations <- matrix(data = NA, nrow = loc_count, ncol = 4)
colnames(target_locations) <- c("gene", "chr", "start", "end")
idx <- 1
#
for (a in 1:loc_count){
  
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
start(gr) <- (start(gr) - 2000)
width(gr) <- (width(gr) + 2000)


## Intersect the targets GRanges with the binding sites list
wg_sites <- mergedMotifs[["sites"]]
#
intersection <- intersect(gr, wg_sites)

## Find overlaps
## YOU WILL WANT TO FIND THE OVERLAPS FOR A RANGE OF EXTENDED VALUES PAST THE TARGETS, GRAPH IT

overlaps <- findOverlaps(gr, wg_sites)




