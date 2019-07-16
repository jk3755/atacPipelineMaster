if(!require(GenomicRanges)){
  cat("Installing GenomicRanges", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")}

if(!require(stats4)){
  cat("Installing stats4", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(BiocGenerics)){
  cat("Installing BiocGenerics", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(parallel)){
  cat("Installing parallel", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(Rsamtools)){
  cat("Installing Rsamtools", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(GenomicAlignments)){
  cat("Installing GenomicAlignments", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(genomation)){
  cat("Installing genomation", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(seqLogo)){
  cat("Installing seqLogo", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(ChIPpeakAnno)){
  cat("Installing ChIPpeakAnno", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(ChIPseeker)){
  cat("Installing ChIPseeker", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(clusterProfiler)){
  cat("Installing clusterProfiler", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(BSgenome.Hsapiens.UCSC.hg38)){
  cat("Installing BSgenome.Hsapiens.UCSC.hg38", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
  cat("Installing TxDb.Hsapiens.UCSC.hg38.knownGene", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(org.Hs.eg.db)){
  cat("Installing org.Hs.eg.db", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(Biostrings)){
  cat("Installing Biostrings", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(MotifDb)){
  cat("Installing MotifDb", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(VariantAnnotation)){
  cat("Installing VariantAnnotation", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

if(!require(ReactomePA)){
  cat("Installing ReactomePA", "\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("stats4")
}

if(!require(ggplot2)){
  cat("Installing ggplot2", "\n")
  install.packages("ggplot2", repos='http://cran.us.r-project.org')
}

if(!require(ggpubr)){
  cat("Installing ggpubr", "\n")
  install.packages("ggpubr", repos='http://cran.us.r-project.org')
}

##
cat("Finished installing packages, touching output file", "\n")
outputPath <- snakemake@output[[1]]
file.create(outputPath)