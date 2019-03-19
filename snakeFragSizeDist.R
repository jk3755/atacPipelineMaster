cat("Loading libraries...", "\n")
#BiocManager::install("ATACseqQC")
suppressMessages(library(ATACseqQC))

##
cat("Setting snakemake vars...", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outputDir <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["mergedsample"]]



##
bamFile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\h508\\wt02a\\H508A-WT-02-REP3.u.bam"
baiFile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\h508\\wt02a\\H508A-WT-02-REP3.u.bai"
##
svg_path <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\h508\\wt02a\\H508A-WT-02-REP3.u.fragsizes.svg"
svg(file = svg_path) # set the filepath for saving the svg figure
##
fragSizeDist(
  bamFiles = bamFile,
  bamFiles.labels = "H508A-WT-02-REP3",
  index = baiFile
)

dev.off()