
## Load libraries
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
footprintData <- aggregateData[["aggregateFootprintData"]]
footprintMetrics <- aggregateData[["aggregateFootprintMetricsData"]]
numGenes <- length(footprintData[,1])

#### ADD RNA EXPRESSION VALUES TO THE DF ####
## Load the RNA expression data from ccle
load("~/git/atacPipelineMaster/atacVIPER/expressionData/ccle/cclecounts.rda")
coadCounts <- cclecounts[["large_intestine_bat1"]]

## Set expression for current cell line
snu61Exp <- coadCounts[,"SNU-61"]
h508Exp <- coadCounts[,"NCI-H508"]
mdst8Exp <- coadCounts[,"MDST8"]
ls1034Exp <- coadCounts[,"LS1034"]

##
curExp <- h508Exp
geneList <- names(curExp)
mappings <- queryMany(geneList, scopes="entrezgene", fields="symbol", species="human")
genemaps <- mappings@listData[["symbol"]]

## Find index of RNA expression values for TFs in this dataset
geneIdx <- which(genemaps %in% footprintData[,"Gene"])

## Add the RNA expression value of each TF to the dataframe
addExp <- c()
for (a in 1:numGenes){
  ## Get the mapping index corresponding to current gene
  curMap <- which(genemaps %in% footprintData[a,1])
  ## If it doesnt exist, set RNAexp to 0
  if (length(curMap) == 0){
    footprintData[a,52] <- 0
  } else {
    footprintData[a,52] <- curExp[[curMap]]
  }} # end for (a in 1:numGenes)


#### ADD VIPER NES SCORE VALUES TO THE DF ####
load("~/git/atacPipelineMaster/atacVIPER/regulons/coad-tcga-regulon.rda")
curVIPER <- msviper(curExp, regul)
curNES <- curVIPER[["es"]][["nes"]]
viperNames <- names(curNES)
viperMap <- queryMany(viperNames, scopes="entrezgene", fields="symbol", species="human")
viperSym <- viperMap@listData[["symbol"]]
viperIdx <- which(viperSym %in% footprintData[,"Gene"])

## Add the VIEPR NES value of each TF to the dataframe
addExp <- c()
for (b in 1:numGenes){
  ## Get the mapping index corresponding to current gene
  curMap <- which(viperSym %in% footprintData[b,1])
  ## If it doesnt exist, set RNAexp to 0
  if (length(curMap) == 0){
    footprintData[b,53] <- 0
  } else {
    footprintData[b,53] <- curNES[[curMap]]
  }} # end for (b in 1:numGenes)


## Transfer to new
H508_data <- footprintData
save(H508_data, file = "C:\\Users\\jsk33\\Desktop\\H508_data.Rdata")

## Set cell line name
cellName <- "MDST8"
footprintData <- MDST8_data

#### GENERATE PLOTS #########################################################################################

## All sites, RNA exp
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("All sites, Depth vs Flank, RNA exp ", cellName))


## All sites, VIPER
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("All sites, Depth vs Flank, VIPER NES ", cellName))


## Bound ratio vs Exp
ggplot(footprintData, aes(boundRatio, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Percent bound sites, RNA exp ", cellName))


## Bound ratio vs VIPER
ggplot(footprintData, aes(boundRatio, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Percent bound sites, VIPER NES ", cellName))


## RNA exp vs Depth
ggplot(footprintData, aes(peak.log2Depth, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("RNA exp vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## RNA exp vs Flank
ggplot(footprintData, aes(peak.log2Flank, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("RNA exp vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## VIPER vs Depth
ggplot(footprintData, aes(peak.log2Depth, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("VIPER vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## VIPER vs Flank
ggplot(footprintData, aes(peak.log2Flank, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("VIPER vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


#### PROMOTERS AND DISTAL ####

## Promoter sites, RNA exp
ggplot(footprintData, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Promoter sites, Depth vs Flank, RNA exp ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band

## Promoter sites, VIPER
ggplot(footprintData, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Promoter sites, Depth vs Flank, VIPER NES ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Promoter sites, RNA exp vs Depth
ggplot(footprintData, aes(promoterPeak.log2Depth, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Promoter sites, RNA exp vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Promoter sites, RNA exp vs Flank
ggplot(footprintData, aes(promoterPeak.log2Flank, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Promoter sites, RNA exp vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Promoter sites, VIPER vs Depth
ggplot(footprintData, aes(promoterPeak.log2Depth, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Promoter sites, VIPER NES vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Promoter sites, VIPER vs Flank
ggplot(footprintData, aes(promoterPeak.log2Flank, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Promoter sites, VIPER NES vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal
## Distal sites, RNA exp
ggplot(footprintData, aes(distalPeak.log2Depth, distalPeak.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Distal sites, Depth vs Flank, RNA exp ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal sites, VIPER
ggplot(footprintData, aes(distalPeak.log2Depth, distalPeak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Distal sites, Depth vs Flank, VIPER NES ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal sites, RNA exp vs Depth
ggplot(footprintData, aes(distalPeak.log2Depth, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Distal sites, RNA exp vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal sites, RNA exp vs Flank
ggplot(footprintData, aes(distalPeak.log2Flank, log2(RNAexp))) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Distal sites, RNA exp vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal sites, VIPER vs Depth
ggplot(footprintData, aes(distalPeak.log2Depth, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Distal sites, VIPER NES vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Distal sites, VIPER vs Flank
ggplot(footprintData, aes(distalPeak.log2Flank, viperNES)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Distal sites, VIPER NES vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Bound promoter
## Bound Promoter sites, RNA exp
ggplot(footprintData, aes(promoterBound.log2Depth, promoterBound.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Bound Promoter sites, Depth vs Flank, RNA exp ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Bound Promoter sites, VIPER
ggplot(footprintData, aes(promoterBound.log2Depth, promoterBound.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Bound Promoter sites, Depth vs Flank, VIPER NES ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## Bound Distal
## Bound Distal sites, RNA exp
ggplot(footprintData, aes(distalBound.log2Depth, distalBound.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Bound Distal sites, Depth vs Flank, RNA exp ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band

## Bound Distal sites, VIPER
ggplot(footprintData, aes(distalBound.log2Depth, distalBound.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  stat_cor(method = "spearman") +
  ggtitle(paste0("Bound Distal sites, Depth vs Flank, VIPER NES ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


####### PLOT COAD MRs, ALL THREE CELL LINES IN SAME PLOT ########

coadmr <- c("CDX2", "TCF7", "MNX1", "POU5F1B", "ESRRA", "HNF4A", "GMEB2", "HOXA3",
            "OVOL1", "ASCL2", "ZSWIM1", "CBFA2T2", "PAX6", "ADNP", "TAF4", "ZMYND8", "ZNF696")

mrIdx <- which(footprintData[,1] %in% coadmr)

## Create new dataframe
MRdata <- data.frame(matrix(vector(), 0, 10,
                                            dimnames=list(c(), c(
                                              "Gene", "Type", "peak.log2Depth", "peak.log2Flank", "RNAexp", "viperNES",
                                              "promoterPeak.log2Depth", "promoterPeak.log2Flank", "bound.log2Depth", "bound.log2Flank"))),
                                     stringsAsFactors=F)

## Transfer info to new dataframe
numPoints <- length(mrIdx)
cellLine <- "MDST8"
cellName <- "MDST8"

## Xfer data
for (d in 1:numPoints){
  
  Idx <- mrIdx[d]
  
  MRdata[d,1] <- footprintData[Idx,1]
  MRdata[d,2] <- cellLine
  MRdata[d,3] <- footprintData[Idx,"peak.log2Depth"]
  MRdata[d,4] <- footprintData[Idx,"peak.log2Flank"]
  MRdata[d,5] <- footprintData[Idx,"RNAexp"]
  MRdata[d,6] <- footprintData[Idx,"viperNES"]
  MRdata[d,7] <- footprintData[Idx,"promoterPeak.log2Depth"]
  MRdata[d,8] <- footprintData[Idx,"promoterPeak.log2Flank"]
  MRdata[d,9] <- footprintData[Idx,"bound.log2Depth"]
  MRdata[d,10] <- footprintData[Idx,"bound.log2Flank"]
}


#### COAD MR PLOTS ##########################################################################################


### Plot
## All sites, RNA exp
ggplot(MRdata, aes(peak.log2Depth, peak.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("COAD MR, Depth vs Flank, RNA exp ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)

ggplot(MRdata, aes(peak.log2Depth, peak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("COAD MR, Depth vs Flank, VIPER NES ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)

#### PROMOTERS
ggplot(MRdata, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Promoters COAD MR, Depth vs Flank, RNA exp ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)

ggplot(MRdata, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Promoters COAD MR, Depth vs Flank, VIPER NES ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)


### Bound sites
ggplot(MRdata, aes(bound.log2Depth, bound.log2Flank, color = log2(RNAexp))) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Bound COAD MR, Depth vs Flank, RNA exp ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)

ggplot(MRdata, aes(bound.log2Depth, bound.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("Bound COAD MR, Depth vs Flank, VIPER NES ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA)
