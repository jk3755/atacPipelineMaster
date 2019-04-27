
## Load libraries and setup
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(mygene))
suppressMessages(library(viper))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(aracne.networks))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
##
options(warn = -1)
options(scipen = 999)

#### GET THE DATA ####
load("C://Users//jsk33//Desktop//aggregate data//H508A-WT-02.aggregated.Rdata") # for testing only
footprintData <- aggregateData[["aggregateFootprintData"]]
footprintMetrics <- aggregateData[["aggregateFootprintMetricsData"]]
numGenes <- length(footprintData[,1])
## Load the regulon
load("~/git/atacPipelineMaster/atacVIPER/regulons/coad-tcga-regulon.rda")
coadInteractome <- regul
## Get symbols of regulator
coadRegulators <- names(coadInteractome)
coadMappings <- queryMany(coadRegulators, scopes="entrezgene", fields="symbol", species="human")
coadMaps <- coadMappings@listData[["symbol"]]
## Which of our TFs with footprints are in the interactome?
matchedRegulatorIdx <- which(coadMaps %in% names(footprintMetrics))
numMatched <- length(matchedRegulatorIdx)

#### CODE TESTING - PROMOTER/DISTAL groups ####
## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
## This modified what metadata is returned from the Txdb database regarding the promoters
#cols <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
cols <- c("gene_id")
promoters <- promoters(txdb, upstream = 1000, downstream = 500, use.names = TRUE, columns = cols)
## Trim the GRanges object to keep standard entries only
scope <- paste0("chr", c(1:22, "X", "Y"))
promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")
#### "unlisting the returned geneids" 
## might need to look at this more closely at some point
## use this code to get a flattened list of the Entrez IDs
#za <- as.matrix(a@listData[["gene_id"]])
promoterIDs <- drop(promoters@elementMetadata@listData[["gene_id"]])


#### Setup a loop to iterate through all regulators also in our ATAC dataset
overlapData <- matrix(data = NA, ncol = 5, nrow = numMatched)
colnames(overlapData) <- c("Total targets", "Peak overlaps", "Bound site overlaps", "Unbound site overlaps", "Percent total overlap")

for (a in 1:numMatched){
  
  tryCatch({
    
    ## Get the name of the current matched regulator
    curName <- coadMaps[matchedRegulatorIdx[a]]
    
    ## Translate that to Entrez ID
    curNameIdx <- which(coadMaps == curName)
    curEntrezID <- coadRegulators[curNameIdx]
    
    ## Pull the targets of the current gene using the Entez ID
    ## For each regulator that has a cooresponding TF footprint,
    ## calculate what percent of the targets are captured by the motif binding sites
    ## Will need to grab the promoter region - extend TSS from -1000/+500 for each gene
    com <- paste0("curTargets <- names(coadInteractome[['", curEntrezID, "']][['tfmode']])")
    eval(parse(text = com))
    
    ## Get the Granges from the peaks of TF footprints for current gene
    com <- paste0("curFootprintsPeaks <- footprintMetrics[['", curName, "']][['peakSites']]")
    eval(parse(text = com))
    ## Get bound ranges
    com <- paste0("curFootprintsBound <- footprintMetrics[['", curName, "']][['boundSites']]")
    eval(parse(text = com))
    ## Get unbound ranges
    com <- paste0("curFootprintsUnbound <- footprintMetrics[['", curName, "']][['unboundSites']]")
    eval(parse(text = com))
    
    ## Which of the returned promoter ranges match to the current target genes? 
    targetPromoterIdx <- which(promoterIDs %in% curTargets)
    ## Subset that Granges for target promoters
    targetPromoters <- promoters[targetPromoterIdx]
    
    ## Find which of the ATAC peak motifs overlap with promoters of any targets for current gene
    targetPeakOverlaps <- findOverlaps(targetPromoters, curFootprintsPeaks, ignore.strand = TRUE)
    targetBoundOverlaps <- findOverlaps(targetPromoters, curFootprintsBound, ignore.strand = TRUE)
    targetUnboundOverlaps <- findOverlaps(targetPromoters, curFootprintsUnbound, ignore.strand = TRUE)
    
    ## Record the extent of the overlap
    overlapData[a,1] <- length(targetPromoters)
    overlapData[a,2] <- length(unique(targetPeakOverlaps@to))
    overlapData[a,3] <- length(unique(targetBoundOverlaps@to))
    overlapData[a,4] <- length(unique(targetUnboundOverlaps@to))
    overlapData[a,5] <- (overlapData[a,2] / overlapData[a,1])*100
    
  }, # end try
  error=function(cond){
    message(cond)
    return(NA)
  },
  finally={})
}

## Plot the overlap data
#plot(density(overlapData[,3]))
#plot(hist(overlapData[,3]))

df <- data.frame(overlapData[,5])
#stip plot
p <- ggplot(df, aes(x=1, y=df[,1])) + geom_jitter(position=position_jitter(0.3)) + theme_bw()
p

p <- ggplot(df, aes(x=1, y=df[,1])) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic()
p

## Modify the COAD interactome according to Footprint data
data("reguloncoad")
write.regulon(regulonblca,file="C:\Users\jsk33\Desktop\aggregate data")


#### MODIFY THE COAD REGULON ####
#### THE FOLLOWING CODE IS ADAPTED FROM CODE ORIGINALLY WRITTEN BY AJAY NAIR ####


# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")

load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
geneId <- "4089" #smad4
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1]==geneId),]#getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1]==geneId)),]# removing the req regulon from interactome
colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
write.table(interactome_tinker, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)
write.table(regulon_req, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
regul <- aracne2regulon("network_coad_tfCotf_3col_tinkered.txt", dset, format = "3col", verbose = TRUE)#regulon
print(regul)
coadViper <- viper(dset, regul, method = "scale") #the viper algorithm
dim(coadViper)


#### THE ABOVE CODE IS ADAPTED FROM CODE ORIGINALLY WRITTEN BY AJAY NAIR ####