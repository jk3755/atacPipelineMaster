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
cat("Loading data", "\n")
load(inputfile)

## USE THIS STRUCTURE FOR HEATMAPS
## first, combine the signals from plus and minus strand
sigs <- parsedSitesInfo[["bfPassPeakSignals"]]
sig_plus <- sigs["+"]
sig_minus <- sigs["-"]
numsites <- length(sig_plus[["+"]][,1])
numbp <- length(sig_plus[["+"]][1,])
## Annotate the combined sublist name which will become the tital of the heatmap plot
plottitle <- paste0(samplename, "_", genename, "_motif", motifnum, "_numsites", parsedSitesInfo[["numbfPassPeakSites"]])
combined <- list()
com <- paste0("combined$", plottitle, " <- matrix(data = NA, nrow = numsites, ncol = numbp)")
eval(parse(text = com))
## combine the signals
com <- paste0("for (i in 1:numsites){combined$", plottitle, "[i,] <- sigs[['+']][i,] + sigs[['-']][i,]}")
eval(parse(text = com))
## set the max value in plot = to max value of combined signals
com <- paste0("maxsig <- max(combined[['", plottitle, "']])")
eval(parse(text = com))
##
svg(file = outputfile) # set the filepath for saving the svg figure
cat("Saving svg footprint image at path:", outputfile, "\n")
ChIPpeakAnno::featureAlignedHeatmap(combined, 
                                    feature.gr=reCenterPeaks(sites,width=num_bp), 
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=num_bp,
                                    upper.extreme = maxsig, # set this to control the heatmap scale
                                    margin = c(0.1, 0.01, 0.05, 0.09),
                                    color=colorRampPalette(c("blue", "white", "red"))(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()
cat("Finished.",)



