
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
    footprintData[a,21] <- 1
  } else {
    footprintData[a,21] <- curExp[[curMap]]
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
    footprintData[b,22] <- 0
  } else {
    footprintData[b,22] <- curNES[[curMap]]
  }} # end for (b in 1:numGenes)

#### GENERATE PLOTS ####

## plot the graph with RNA expression
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = RNAexp)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

## plot the graph with viper values
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank, color = viperNES)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")

## Plot RNA exp and NES against flank and depth individually



## Plot 1 - Depth vs Flank

# Peaks
ggplot(footprintData, aes(peak.log2Depth, peak.log2Flank)) + 
  geom_point()

# Bound
ggplot(footprintData, aes(bound.log2Depth, bound.log2Flank)) + 
  geom_point()

# Unbound
ggplot(footprintData, aes(unbound.log2Depth, unbound.log2Flank)) + 
  geom_point()

## Plot 2 - Num unboud vs bound sites
ggplot(footprintData, aes(log(numBoundSites), log(numUnboundSites))) + 
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) # se controls the confidence band



#### PLOT SINGLE TF SITES ####
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