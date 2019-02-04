## Load libraries
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(rtracklayer))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(MotifDb))

## Set snakemake variables
cat("Setting snakemake vars...", "\n")
inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["mergedsample"]]
genename <- snakemake@wildcards[["gene"]]
nummotif <- snakemake@wildcards[["nummotif"]]
dirpath <- snakemake@wildcards[["path"]]

##
widths <- c()
totalsites <- 0
totalBp <- 0
signalsMerged <- list()

## Load the individual motifs, do one by one for memory considerations
cat("Loading data...", "\n")
for (a in 1:nummotif){
  
  filepath <- gsub("parsed.done.txt", paste0("motif", a, ".info.Rdata"), inputfile)
  load(filepath)
  # signals
  com <- paste0("signals", a, " <- parsedSitesInfo[['bfPassPeakSignals']]")
  eval(parse(text = com))
  # sites
  com <- paste0("sites", a, " <- parsedSitesInfo[['bfPassPeakSites']]")
  eval(parse(text = com))
  # number of sites
  com <- paste0("numsites", a, " <- length(signals", a, "[['+']][,1])")
  eval(parse(text = com))
  # motif width
  com <- paste0("motifwidth", a, " <- sites", a, "@ranges@width[1]")
  eval(parse(text = com))
  # add
  com <- paste0("widths <- c(widths, motifwidth", a, ")")
  eval(parse(text = com))
  com <- paste0("totalsites <- (totalsites + length(signals", a, "[['+']][,1]))")
  eval(parse(text = com))
  
}

## Must set the width to the length of the smallest vector and set totalbp
minwidth <- min(widths)
totalbp <- (minwidth+200)


## To center the Granges, change the width of all sites to the size of smallest motif width
cat("Centering GRanges...", "\n")
for (b in 1:nummotif){
  com <- paste0("width(sites", b, ") <- minwidth")
  eval(parse(text = com))
}

## Now that all granges are centered together, can simply remove trailing values
cat("Combining and merging signals...", "\n")
for (c in 1:nummotif){
  
  com <- paste0("sigplus <- signals", c, "[['+']]")
  eval(parse(text = com))
  com <- paste0("sigminus <- signals", c, "[['-']]")
  eval(parse(text = com))
  #
  com <- paste0("sigs", c, "new <- matrix(data = NA, nrow = length(sigplus[,1]), ncol = totalbp)")
  eval(parse(text = com))
  #
  com <- paste0("for (d in 1:length(sigplus[,1])){for (e in 1:totalbp){sigs", c, "new[c,d] <- (sigplus[c,d] + sigminus[c,d])")
  eval(parse(text = com))
}

##
mergePlus <- matrix(data = NA, nrow = totalsites, ncol = totalbp)
##
mergeSignals <- sigs1new
for (e in 2:nummotif){
  
  com <- paste0("mergeSignals <- rbind(mergeSignals, sigs", e, "new)")
  eval(parse(text = com))

}

##
sigs <- list()
sigs$signal <- mergeSignals

## FIX ME
cat("Merging Granges...", "\n")
mergesites <- c(sites1,sites2,sites3,sites4)


##
cat("Generating heatmap...", "\n")
heatmappath <- paste0("dirpath/merged_motifs/", samplename, ".", genename, ".merged.heatmap.svg")
svg(file = heatmappath)
ChIPpeakAnno::featureAlignedHeatmap(sigs,
                                    feature.gr=reCenterPeaks(mergesites,width=totalbp), 
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=totalbp,
                                    #upper.extreme = maxsig, # set this to control the heatmap scale
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("blue", "white", "yellow", "red"), bias=3)(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()

##
cat("Saving merged info...", "\n")
merged <- list()
merged$sigs <- sigs
merged$sites <- mergedsites
save(merged, file = outputfile)



