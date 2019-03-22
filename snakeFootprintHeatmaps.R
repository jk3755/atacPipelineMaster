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
inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
genename <- snakemake@wildcards[["gene"]]
motifnum <- snakemake@wildcards[["motif"]]

## Loading data
cat("Loading data...", "\n")
load(inputfile)

## USE THIS STRUCTURE FOR HEATMAPS
## first, combine the signals from plus and minus strand
sigs <- parsedSitesInfo[["bfPassPeakSignals"]]
sig_plus <- sigs["+"]
sig_minus <- sigs["-"]
numsites <- length(sig_plus[["+"]][,1])
numbp <- length(sig_plus[["+"]][1,])
sites <- parsedSitesInfo[["bfPassPeakSites"]]

##
temp <- matrix(data = NA, nrow = numsites, ncol = numbp)
for (i in 1:numsites){temp[i,] <- sigs[['+']][i,] + sigs[['-']][i,]}

## scale each row individually
for (a in 1:numsites){
  maxsig <- max(temp[a,])
  for (b in 1:numbp){temp[a,b] <- (temp[a,b] / maxsig)}}
maxsig <- 1

## invert signals
for (a in 1:numsites){for (b in 1:numbp){temp[a,b] <- (1-temp[a,b])}}


## Annotate the combined sublist name which will become the tital of the heatmap plot
plottitle <- paste0(genename, "_motif", motifnum, "_numsites", parsedSitesInfo[["numbfPassPeakSites"]])
combined <- list()
com <- paste0("combined$", plottitle, " <- temp")
eval(parse(text = com))


##
svg(file = outputfile) # set the filepath for saving the svg figure
cat("Saving svg footprint image at path:", outputfile, "\n")

## Margin controls
# margin(a,b,c,d)
# a = size of graph from top to bottom, higher value = smaller. default = 0.1
# b = size of graph from left to right, higher value = smaller. default = 0.005
# c = flips x axis?
# d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
# good settings for ATACseq = c(0.1,0.005,0.05,0.2)
# bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values

##
ChIPpeakAnno::featureAlignedHeatmap(combined,
                                    feature.gr=reCenterPeaks(sites,width=numbp),
                                    upper.extreme = maxsig, # set this to control the heatmap scale
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=numbp,
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias=0.9)(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()
cat("Finished.", "\n")



