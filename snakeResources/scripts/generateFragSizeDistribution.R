## Set snakemake variables
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
outPath <- snakemake@output[[1]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
svgOut1 <- gsub("operations/metrics", "metrics", outPath)
svgOut1 <- gsub("fragsizes.done", "fragsize1.svg", svgOut1)
svgOut2 <- gsub("fragsize1", "fragsize2", svgOut1)
svgOut3 <- gsub("fragsize1", "fragsize3", svgOut1)
svgOut4 <- gsub("fragsize1", "fragsize4", svgOut1)

## Report
cat("Bam filepath:", bamPath, "\n")
cat("Bai filepath:", baiPath, "\n")
cat("Output filepath:", outPath, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")
cat("Output for svg plot 1:", svgOut1, "\n")
cat("Output for svg plot 2:", svgOut2, "\n")
cat("Output for svg plot 3:", svgOut3, "\n")
cat("Output for svg plot 4:", svgOut4, "\n")

## Load libraries
cat("Loading libraries", "\n")
suppressMessages(library(GenomicAlignments))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

## Load the reads
cat("Loading reads", "\n")
atacReads <- readGAlignmentPairs(bamPath, param = ScanBamParam(mapqFilter = 1, 
                                                                 flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
                                                                 "mapq", "isize")))
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

## Plot 1
cat("Generating fragment size distribution 1", "\n")
svg(file = svgOut1)
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                            Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
                                                            geom_line()
fragLenPlot + theme_bw()
dev.off()

## Plot 2
cat("Generating fragment size distribution 2", "\n")
svg(file = svgOut2)
fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
dev.off()

## Plot 3
cat("Generating fragment size distribution 3", "\n")
svg(file = svgOut3)
fragLenPlot + geom_vline(xintercept = c(180, 247),
                         colour = "red") + geom_vline(xintercept = c(315, 437),
                         colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
dev.off()

## Plot 4
cat("Generating fragment size distribution 4", "\n")
svg(file = svgOut4)
fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 247),
                                colour = "red") + geom_vline(xintercept = c(315, 437),
                                colour = "darkblue") +
                                geom_vline(xintercept = c(100),
                                colour = "darkgreen") + theme_bw()
dev.off()

##
cat("Finished, touching snakemake flag file", "\n")
file.create(outPath)