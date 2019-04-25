
## Load libraries
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(mygene))
suppressMessages(library(viper))
suppressMessages(library(ggrepel))

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
curExp <- mdst8Exp
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
MDST8_data <- footprintData
save(MDST8_data, file = "C:\\Users\\jsk33\\Desktop\\MDST8_data.Rdata")


## Set cell line name
cellName <- "MDST8"




##### TO ADD A NEW COLUMN TO DF ####
## bound ratio
#x <- footprintData[,4]/footprintData[,3]
#footprintData["Ratio"] <- x



#### COAD MR PLOTS ##########################################################################################


## All sites, RNA exp
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle(paste0("COAD MR Depth vs Flank, RNA exp ", cellName)) +
  geom_label_repel(aes(label = Gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   label.size = NA) +
  theme_classic()

## All sites, VIPER
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("COAD MR Depth vs Flank, VIPER NES, MDST8")

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
  ggtitle(paste0("Percent bound sites, RNA exp ", cellName))


## Bound ratio vs VIPER
ggplot(footprintData, aes(boundRatio, viperNES)) + 
  geom_point() +
  ggtitle(paste0("Percent bound sites, VIPER NES ", cellName))


## RNA exp vs Depth
ggplot(footprintData, aes(peak.log2Depth, log2(RNAexp))) + 
  geom_point() +
  ggtitle(paste0("RNA exp vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## RNA exp vs Flank
ggplot(footprintData, aes(peak.log2Flank, log2(RNAexp))) + 
  geom_point() +
  ggtitle(paste0("RNA exp vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## VIPER vs Depth
ggplot(footprintData, aes(peak.log2Depth, viperNES)) + 
  geom_point() +
  ggtitle(paste0("VIPER vs Depth ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


## VIPER vs Flank
ggplot(footprintData, aes(peak.log2Flank, viperNES)) + 
  geom_point() +
  ggtitle(paste0("VIPER vs Flank ", cellName)) +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band


#### PROMOTERS AND DISTAL ####

## Promoter sites, RNA exp
ggplot(footprintData, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Promoter sites, Depth vs Flank, RNA exp, MDST8")

## Promoter sites, VIPER
ggplot(footprintData, aes(promoterPeak.log2Depth, promoterPeak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Promoter sites, Depth vs Flank, VIPER NES, MDST8")


## Promoter sites, RNA exp vs Depth
ggplot(footprintData, aes(promoterPeak.log2Depth, RNAexp)) + 
  geom_point() +
  ggtitle("Promoter sites, RNA exp vs Depth, MDST8")


## Promoter sites, RNA exp vs Flank
ggplot(footprintData, aes(promoterPeak.log2Flank, RNAexp)) + 
  geom_point() +
  ggtitle("Promoter sites, RNA exp vs Flank, MDST8")


## Promoter sites, VIPER vs Depth
ggplot(footprintData, aes(promoterPeak.log2Depth, viperNES)) + 
  geom_point() +
  ggtitle("Promoter sites, VIPER NES vs Depth, MDST8")


## Promoter sites, VIPER vs Flank
ggplot(footprintData, aes(promoterPeak.log2Flank, viperNES)) + 
  geom_point() +
  ggtitle("Promoter sites, VIPER NES vs Flank, MDST8")


## Distal
## Distal sites, RNA exp
ggplot(footprintData, aes(distalPeak.log2Depth, distalPeak.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Distal sites, Depth vs Flank, RNA exp, MDST8")

## Distal sites, VIPER
ggplot(footprintData, aes(distalPeak.log2Depth, distalPeak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Distal sites, Depth vs Flank, VIPER NES, MDST8")


## Distal sites, RNA exp vs Depth
ggplot(footprintData, aes(distalPeak.log2Depth, RNAexp)) + 
  geom_point() +
  ggtitle("Distal sites, RNA exp vs Depth, MDST8")


## Distal sites, RNA exp vs Flank
ggplot(footprintData, aes(distalPeak.log2Flank, RNAexp)) + 
  geom_point() +
  ggtitle("Distal sites, RNA exp vs Flank, MDST8")


## Distal sites, VIPER vs Depth
ggplot(footprintData, aes(distalPeak.log2Depth, viperNES)) + 
  geom_point() +
  ggtitle("Distal sites, VIPER NES vs Depth, MDST8")


## Distal sites, VIPER vs Flank
ggplot(footprintData, aes(distalPeak.log2Flank, viperNES)) + 
  geom_point() +
  ggtitle("Distal sites, VIPER NES vs Flank, MDST8")


## Bound promoter
## Bound Promoter sites, RNA exp
ggplot(footprintData, aes(promoterBound.log2Depth, promoterBound.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Bound Promoter sites, Depth vs Flank, RNA exp, MDST8")

## Bound Promoter sites, VIPER
ggplot(footprintData, aes(promoterBound.log2Depth, promoterBound.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Bound Promoter sites, Depth vs Flank, VIPER NES, MDST8")


## Bound Distal
## Bound Distal sites, RNA exp
ggplot(footprintData, aes(distalBound.log2Depth, distalBound.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Bound Distal sites, Depth vs Flank, RNA exp, MDST8")

## Bound Distal sites, VIPER
ggplot(footprintData, aes(distalBound.log2Depth, distalBound.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Bound Distal sites, Depth vs Flank, VIPER NES, MDST8")


####### PLOT COAD MRs, ALL THREE CELL LINES IN SAME PLOT ########

## Create new dataframe
aggregateFootprintData <- data.frame(matrix(vector(), 0, 4,
                                            dimnames=list(c(), c(
                                              "Gene", "Type", "peak.log2Depth", "peak.log2Flank", "RNAexp", "viperNES"))),
                                     stringsAsFactors=F)



## Transfer info to new dataframe
numPoints <- length(footprintData[,1])
cellLine <- "MDST8"

## Xfer data
for (x in 1:numPoints){
  
  aggregateFootprintData[x,1] <- footprintData[x,1]
  aggregateFootprintData[x,2] <- cellLine
  aggregateFootprintData[x,3] <- footprintData["peak.log2Depth"]
  aggregateFootprintData[x,4] <- footprintData["peak.log2Flank"]
  aggregateFootprintData[x,5] <- footprintData["RNAexp"]
  aggregateFootprintData[x,6] <- footprintData["viperNES"]
  
}


geom_point(aes(colour = factor(cyl)))

## All cell lines, RNA exp
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Depth vs Flank, RNA exp, MDST8")


## All cell lines, VIPER
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red") +
  ggtitle("Depth vs Flank, VIPER NES, MDST8")
















ggplot(footprintData, aes(Ratio, viperNES)) + 
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band

#### PLOT SINGLE TF SITES ####
ggplot(footprintData, aes(log(numBoundSites), log(numUnboundSites))) + 
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band

ggplot(footprintData, aes(log(numBoundSites), log(numUnboundSites))) + 
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band

peakMetrics <- data.frame(footprintMetrics[["ATF1"]][["rawPeakFootprintMetrics"]])
boundMetrics <- data.frame(footprintMetrics[["ATF1"]][["boundSitesMetrics"]])
unboundMetrics <- data.frame(footprintMetrics[["ATF1"]][["unboundSitesMetrics"]])

# Flank vs motif signal
ggplot(peakMetrics, aes(Flanking, Motif)) + 
  geom_point() +
  ggtitle("Peaks")

ggplot(boundMetrics, aes(Flanking, Motif)) + 
  geom_point() +
  ggtitle("Bound Sites")

ggplot(unboundMetrics, aes(Flanking, Motif)) + 
  geom_point() +
  ggtitle("Unbound Sites")


### Flanking acc vs depth
ggplot(peakMetrics, aes(Flanking.Accessibility, Footprint.Depth)) + 
  geom_point() +
  ggtitle("Peaks, Flanking vs. Depth")

### Flanking acc vs depth
ggplot(boundMetrics, aes(Flanking.Accessibility, Footprint.Depth)) + 
  geom_point() +
  ggtitle("Bound, Flanking vs. Depth")

### Flanking acc vs depth
ggplot(unboundMetrics, aes(Flanking.Accessibility, Footprint.Depth)) + 
  geom_point() +
  ggtitle("Unbound, Flanking vs. Depth")

### SHOW BOTH POPULATIONS IN ONE GRAPH, COLOR POINTS BY BOUND/UNBOUND

#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(aggregateFootprintData, aes(peak.log2depth, peak.log2flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")


#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")


#### Generate the plots ####
## Plot 1
## Footprint depth against flanking accessibility
## Data points colored by RNA expression level
ggplot(dfFootprints, aes(depth, flank, color=exp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

## plot the graph with viper values
ggplot(dfFootprints, aes(depth, flank, color=nex)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")



## plot footprint depth against viper activity
ggplot(dfFootprints, aes(nex, depth, color=flank)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")