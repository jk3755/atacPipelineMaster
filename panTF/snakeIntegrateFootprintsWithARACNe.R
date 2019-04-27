
## Load libraries and setup
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(mygene))
suppressMessages(library(viper))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
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
coad_interactome <- regul

## Get symbols of regulator
coadRegulators <- names(coad_interactome)
coadMappings <- queryMany(coadRegulators, scopes="entrezgene", fields="symbol", species="human")
coadMaps <- coadMappings@listData[["symbol"]]


## Which of our TFs with footprints are in the interactome?
matchedRegulatorIdx <- which(coadMaps %in% names(footprintMetrics))
numMatched <- length(matchedRegulatorIdx)


#### Setup a loop to iterate through all regulators also in our ATAC dataset
overlapData <- matrix(data = NA, ncol = 2, nrow = numMatched)

for (a in 1:numMatched){
  
  ## Get the name of the current matched regulator
  curName <- names(footprintMetrics)[matchedRegulatorIdx[1]]
  
  ## Translate that to Entrez ID
  curNameIdx <- which(coadMaps == curName)
  curEntrezID <- coadRegulators[curNameIdx]
  
  ## Pull the targets of the current gene using the Entez ID
  ## For each regulator that has a cooresponding TF footprint,
  ## calculate what percent of the targets are captured by the motif binding sites
  ## Will need to grab the promoter region - extend TSS from -1000/+500 for each gene
  com <- paste0("curTargets <- names(coad_interactome[['", curEntrezID, "']][['tfmode']])")
  eval(parse(text = com))
  
  
}

a <- which(names(footprintMetrics) == coadMaps[2])





#### CODE TESTING - PROMOTER/DISTAL groups ####
## Pull promoters from txdb, define promoter region as -1000/+100 in accordance with TCGA paper
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoters <- promoters(txdb, upstream = 1000, downstream = 500, use.names = TRUE, columns = cols)
## Trim the GRanges object to keep standard entries only
scope <- paste0("chr", c(1:22, "X", "Y"))
promoters <- keepStandardChromosomes(promoters, pruning.mode="coarse")
promoters <- keepSeqlevels(promoters, scope, pruning.mode="coarse")

## Subset based on the overlaps
promoterOverlaps <- findOverlaps(promoters, peakSites, ignore.strand = TRUE)
promoterIdx <- unique(promoterOverlaps@to)
##
promoterPeakSites <- peakSites[promoterIdx]
distalPeakSites <- peakSites[-promoterIdx]

## This modified what metadata is returned from the Txdb database regarding the promoters
#cols <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
cols <- c("gene_id")


#### "unlisting the returned geneids" 
## might need to look at this more closely at some point
#
ll <- drop(a@listData[["gene_id"]])
#za <- as.matrix(a@listData[["gene_id"]])


##
id <- which(targets %in% ll)
id <- which(ll %in% targets)


##
newPromoters <- promoters[id]

##
peakSites <- footprintMetrics[["ZBTB33"]][["peakSites"]]
ov <- findOverlaps(newPromoters, peakSites, ignore.strand = TRUE)
