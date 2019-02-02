cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(MotifDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(Rsamtools))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(BiocGenerics))
suppressMessages(library(parallel))

##
cat("Setting snakemake vars...", "\n")
bampath <- snakemake@input[[1]]
baipath <- snakemake@input[[2]]
sitespath <- snakemake@input[[3]]
final_outpath <- snakemake@output[[1]]
samplename <- snakemake@wildcards[["sample"]]
genename <- snakemake@wildcards[["gene"]]
current_chr <- snakemake@wildcards[["chr"]]
sample_path <- snakemake@wildcards[["path"]]

##
cat("Loading binding sites and PWM...", "\n")
load(sitespath)
num_motif <- length(bindingSites)