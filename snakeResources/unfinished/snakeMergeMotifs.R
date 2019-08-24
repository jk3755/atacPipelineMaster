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
totalbp <- 0
mergeSignals <- list()

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
  #
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
cat("Combining plus and minus strand signals and annotating by total row signal...", "\n")
for (c in 1:nummotif){
  
  ##
  com <- paste0("sigplus <- signals", c, "[['+']]")
  eval(parse(text = com))
  
  ##
  com <- paste0("sigminus <- signals", c, "[['-']]")
  eval(parse(text = com))
  
  ##
  com <- paste0("sigs", c, "new <- matrix(data = NA, nrow = length(sigplus[,1]), ncol = totalbp)")
  eval(parse(text = com))
  
  ##
  com <- paste0("for (d in 1:length(sigplus[,1])){for (e in 1:totalbp){sigs", c, "new[d,e] <- (sigplus[d,e] + sigminus[d,e])}}")
  eval(parse(text = com))
  
}


## Merge the signals with rbind
cat("Merging signals...", "\n")
mergeSignals <- sigs1new
for (e in 2:nummotif){
  com <- paste0("mergeSignals <- rbind(mergeSignals, sigs", e, "new)")
  eval(parse(text = com))
}

## Find row totals from merged signals
cat("Calculating row totals...", "\n")
rowtotals <- c()
nsites <- length(mergeSignals[,1])
com <- paste0("for (x in 1:nsites){rowtotals[x] <- sum(mergeSignals[x,])}")
eval(parse(text = com))

##
cat("Merging GRanges...", "\n")
temp1 <- paste0("sites", c(1:nummotif), ",")
temp2 <- gsub(paste0(nummotif,","), nummotif, temp1)
temp3 <- paste(temp2, collapse = " ")
com <- paste0("mergedsites <- c(", temp3, ")")
eval(parse(text = com))
## Add annotation column for row total signal
cat("Adding row totals annotation to GRanges...", "\n")
mergedsites@elementMetadata@listData$rowtotals <- rowtotals


## set the max value
cat("Finding max signal...", "\n")
maxsig <- max(mergeSignals)
## normalize all values to max signal
cat("Normalizing signals...", "\n")
com <- paste0("for (f in 1:totalsites){for (g in 1:totalbp){mergeSignals[f,g] <- (mergeSignals[f,g]/maxsig)}}")
eval(parse(text = com))
## reset maxsig for plotting purposes, should always equal 1 now
maxsig2 <- max(mergeSignals)

## This format required for heatmap generation
# generate plot title
plottitle <- paste0(genename, "_mergedMotifs", "_numSites", totalsites)
sigs <- list()
com <- paste0("sigs$", plottitle, " <- mergeSignals")
eval(parse(text = com))

##
## For merged motifs that have a ton of sites, like >100,000, can increase the (100) next to colorRampPallette, increase to 1000
cat("Generating heatmap...", "\n")
heatmappath <- paste0(dirpath, "merged_motifs/", samplename, ".", genename, ".merged.heatmap.svg")
svg(file = heatmappath)
ChIPpeakAnno::featureAlignedHeatmap(sigs,
                                    feature.gr=reCenterPeaks(mergedsites,width=totalbp), 
                                    annoMcols="rowtotals",
                                    sortBy="rowtotals",
                                    n.tile=totalbp,
                                    upper.extreme = maxsig2,
                                    margin = c(0.1, 0.005, 0.05, 0.2),
                                    color=colorRampPalette(c("blue", "white", "yellow", "red"), bias=3)(100),
                                    gp = gpar(fontsize=10),
                                    newpage = TRUE)
dev.off()

##
cat("Processed a total of ", totalsites, "sites from ", nummotif, "unique binding motifs", "\n")
cat("Saving merged info...", "\n")
mergedMotifs <- list()
mergedMotifs$sigs <- sigs
mergedMotifs$sites <- mergedsites
mergedMotifs$gene <- genename
mergedMotifs$totalmotif <- nummotif
mergedMotifs$totalsites <- totalsites
mergedMotifs$totalbp <- totalbp
mergedMotifs$maxsig <- maxsig
save(mergedMotifs, file = outputfile)



