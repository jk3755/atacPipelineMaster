#### Set snakemake variables ####
cat("Setting snakemake variables", "\n")
bamPath <- snakemake@input[[1]]
baiPath <- snakemake@input[[2]]
sampleName <- snakemake@wildcards[["sample"]]
sampleRep <- snakemake@wildcards[["repnum"]]
refGenome <- snakemake@wildcards[["refgenome"]]

#### Output paths ####
svgOut1 <- snakemake@output[[1]]
svgOut2 <- snakemake@output[[2]]
svgOut3 <- snakemake@output[[3]]
svgOut4 <- snakemake@output[[4]]

#### Report ####
cat("Bam filepath:", bamPath, "\n")
cat("Bai filepath:", baiPath, "\n")
cat("Sample name:", sampleName, "\n")
cat("Sample rep:", sampleRep, "\n")
cat("Reference genome:", refGenome, "\n")
cat("Output for svg plot 1:", svgOut1, "\n")
cat("Output for svg plot 2:", svgOut2, "\n")
cat("Output for svg plot 3:", svgOut3, "\n")
cat("Output for svg plot 4:", svgOut4, "\n")

#### Load libraries ####
cat("Loading libraries", "\n")
suppressMessages(library(GenomicAlignments))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

#### Load the reads ####
cat("Loading reads", "\n")
atacReads <- readGAlignmentPairs(
                                bamPath,
                                param = ScanBamParam(mapqFilter = 1,
                                                     flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
                                                     what = c("qname", "mapq", "isize")))

####
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

#### Plot 1 ####
cat("Generating fragment size distribution 1", "\n")
svg(filename = svgOut1)
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                            Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
                                                            geom_line()
fragLenPlot + theme_bw()
dev.off()

#### Plot 2 ####
cat("Generating fragment size distribution 2", "\n")
svg(filename = svgOut2)
fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
dev.off()

#### Plot 3 ####
cat("Generating fragment size distribution 3", "\n")
svg(filename = svgOut3)
fragLenPlot + geom_vline(xintercept = c(180, 247),
                         colour = "red") + geom_vline(xintercept = c(315, 437),
                         colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
dev.off()

#### Plot 4 ####
cat("Generating fragment size distribution 4", "\n")
svg(filename = svgOut4)
fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 247),
                                colour = "red") + geom_vline(xintercept = c(315, 437),
                                colour = "darkblue") +
                                geom_vline(xintercept = c(100),
                                colour = "darkgreen") + theme_bw()
dev.off()
