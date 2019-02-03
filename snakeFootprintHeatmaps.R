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
## Annotate the combined sublist name which will become the tital of the heatmap plot
plottitle <- paste0(genename, "_motif", motifnum, "_numsites", parsedSitesInfo[["numbfPassPeakSites"]])
combined <- list()
com <- paste0("combined$", plottitle, " <- matrix(data = NA, nrow = numsites, ncol = numbp)")
eval(parse(text = com))
## combine the signals
com <- paste0("for (i in 1:numsites){combined$", plottitle, "[i,] <- sigs[['+']][i,] + sigs[['-']][i,]}")
eval(parse(text = com))
## Add annotation column for row total signal
sites <- parsedSitesInfo[["bfPassPeakSites"]]
rowtotals <- c()
for (x in 1:numsites){rowtotals[x] <- sum(combined[[1]][x,])}
sites@elementMetadata@listData$rowtotal <- rowtotals
## set the max value in plot = to max value of combined signals
com <- paste0("maxsig <- max(combined[['", plottitle, "']])")
eval(parse(text = com))
## normalize all values to max signal
com <- paste0("for (a in 1:numsites){for (b in 1:numbp){combined$", plottitle, "[a,b] <- (combined$", plottitle, "[a,b]/maxsig)}}")
eval(parse(text = com))
com <- paste0("maxsig <- max(combined[['", plottitle, "']])")
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
ChIPpeakAnno::featureAlignedHeatmap(combined,
                                    feature.gr=reCenterPeaks(sites,width=numbp), 
                                    annoMcols="rowtotal",
                                    sortBy="rowtotal",
                                    n.tile=numbp,
                                    upper.extreme = maxsig, # set this to control the heatmap scale
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("blue", "white", "red"))(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()
cat("Finished.", "\n")



