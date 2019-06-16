##
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(MotifDb))

##
cat("Setting snakemake vars...", "\n")
sample1Data <- snakemake@input[[1]]
sample2Data <- snakemake@input[[2]]
outputPath <- snakemake@output[[1]]
sampleName1 <- snakemake@wildcards[["xsample1"]]
sampleName2 <- snakemake@wildcards[["xsample2"]]
geneName <- snakemake@wildcards[["gene"]]
motifNum <- snakemake@wildcards[["motif"]]

##
cat("Analyzing samples", sampleName1, sampleName2, "for gene", geneName, "motif number", motifNum, "\n")

## Set PWM for this motif
PWMin <- pwm2pfm(PWM)

## Load in the data from both samples
cat("Loading data for first sample from path:", sample1Data, "\n")
load(sample1Data)
cat("Loading data for second sample from path:", sample2Data, "\n")
load(sample2Data)

## Sample 1 setup
sample1Signals <-
sample1Sites <- 

## Sample 2 setup
sample2Signals <-
sample2Sites <- 

## Find GRanges intersection for bf peak passing sites for both samples
xSitesOverlap <- findOverlaps(sample1Sites, sample2Sites)
#xSites <- GenomicRanges::intersect(sample1Sites, sample2Sites)

## Subset the signals and sites objects for both samples based on the intersection
idx1 <- xSitesOverlap@from
idx2 <- xSitesOverlap@to


## Get common sites
sample1xSites <- sample1Sites[idx1]
sample2xSites <- sample2Sites[idx2]


## Get common signals
sample1xSignals <- sample1Signals[idx1,]
sample2xSignals <- sample2Signals[idx2,]




## Make sure the data is in the same order for both samples



## Generate and save footprint plots
svgPlot1 <- paste0("xsample/", sampleName1, ".", sampleName2, "/footprints/", sampleName1, ".", geneName, ".", motifNum, ".bfpeak.xsites.svg")
cat("Saving footprint plot for sample 1 at path:", svgPlot1, "\n")

##
svg(file = svgPlot1) # set the filepath for saving the svg figure
ATACseqQC:::plotFootprints(prof,
                           Mlen=motifWidth, motif=PWMin)
dev.off()



## Generate and save heatmap footprint plots

cat("Generating heatmap...", "\n")
cat("1", "\n")
heatsigs <- combinedbfPassPeakSignal
cat("2", "\n")
heatnumsites <- numbfPassPeakSites
cat("3", "\n")
heatnumbp <- length(combinedbfPassPeakSignal[1,])
cat("4", "\n")
heatsites <- bfPassPeakSites


#sigs <- parsedSitesInfo[["bfPassPeakSignals"]]
#sig_plus <- sigs["+"]
#sig_minus <- sigs["-"]
#numsites <- length(sig_plus[["+"]][,1])
#numbp <- length(sig_plus[["+"]][1,])
#sites <- parsedSitesInfo[["bfPassPeakSites"]]

##
#temp <- matrix(data = NA, nrow = numsites, ncol = numbp)
#for (i in 1:numsites){temp[i,] <- sigs[['+']][i,] + sigs[['-']][i,]}

## scale each row individually
cat("5", "\n")
for (a in 1:heatnumsites){
  maxsig <- max(heatsigs[a,])
  for (b in 1:heatnumbp){heatsigs[a,b] <- (heatsigs[a,b] / maxsig)}}
maxsig <- 1

## invert signals
cat("6", "\n")
for (a in 1:heatnumsites){for (b in 1:heatnumbp){heatsigs[a,b] <- (1-heatsigs[a,b])}}

## Annotate the combined sublist name which will become the title of the heatmap plot
cat("7", "\n")
plottitle <- paste0(genename, "_motif", x, "_numsites", heatnumsites)
combined <- list()
com <- paste0("combined$", plottitle, " <- heatsigs")
eval(parse(text = com))

##
cat("8", "\n")
heatsvg_path <- paste0(dirpath, "footprints/heatmaps/", samplename, ".", genename, ".", "motif", x, ".bfpeak.sites.heatmap.svg")
svg(file = heatsvg_path) # set the filepath for saving the svg figure
cat("Saving svg footprint image at path:", heatsvg_path, "\n")

## Margin controls
# margin(a,b,c,d)
# a = size of graph from top to bottom, higher value = smaller. default = 0.1
# b = size of graph from left to right, higher value = smaller. default = 0.005
# c = flips x axis?
# d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
# good settings for ATACseq = c(0.1,0.005,0.05,0.2)
# bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values

##
cat("9", "\n")
ChIPpeakAnno::featureAlignedHeatmap(combined,
                                    feature.gr=reCenterPeaks(heatsites,width=heatnumbp),
                                    upper.extreme = maxsig, # set this to control the heatmap scale
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=heatnumbp,
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias=0.9)(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()
############################################################









