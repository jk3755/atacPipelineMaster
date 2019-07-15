#### Install libraries, if necessary ####

## Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges", suppressUpdates = TRUE)
biocLite("stats4", suppressUpdates = TRUE)
biocLite("BiocGenerics", suppressUpdates = TRUE)
biocLite("parallel", suppressUpdates = TRUE)
biocLite("Rsamtools", suppressUpdates = TRUE)
biocLite("GenomicAlignments", suppressUpdates = TRUE)
biocLite("genomation", suppressUpdates = TRUE)
biocLite("seqLogo", suppressUpdates = TRUE)
biocLite("ChIPpeakAnno", suppressUpdates = TRUE)
biocLite("ChIPseeker", suppressUpdates = TRUE)
biocLite("clusterProfiler", suppressUpdates = TRUE)
biocLite("BSgenome.Hsapiens.UCSC.hg38", suppressUpdates = TRUE)
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene", suppressUpdates = TRUE)
biocLite("org.Hs.eg.db", suppressUpdates = TRUE)
biocLite("Biostrings", suppressUpdates = TRUE)
biocLite("MotifDb", suppressUpdates = TRUE)
biocLite("VariantAnnotation", suppressUpdates = TRUE)
biocLite("ReactomePA", suppressUpdates = TRUE)

## Base
install.packages("ggplot2", repos='http://cran.us.r-project.org')
install.packages("ggpubr", repos='http://cran.us.r-project.org')
