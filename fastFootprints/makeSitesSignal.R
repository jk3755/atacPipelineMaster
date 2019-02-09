
## Load libraries
gc()
cat("Loading libraries...", "\n")
suppressMessages(library(ATACseqQC))
suppressMessages(library(MotifDb))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

## Snakemake variables
cat("Loading snakemake variables...", "\n")
insmatrix_path <- snakemake@input[[1]]
sites_path <- snakemake@input[[2]]
out_filepath <- snakemake@output[[1]]
current_chr <- snakemake@wildcards[["chr"]]


## Genomic sites
cat("Loading genomic sites...", "\n")
load("C:\\Users\\jsk33\\Documents\\atac\\coad_mr_sites\\CDX2.sites.Rdata")
load(sites_path)
sites <- data.frame(bindingSites[[1]][["sites"]])
rm(bindingSites, PWMlist, temp, b, num_pwm, outpath, genome)
motif_size <- length(PWM[1,])
# subset for current xsome
sites <- sites[which(sites[,1] == current_chr),]
# number of sites under consideration
num_sites <- length(sites[,1])
# declare the matrix
sites_signals <- matrix(data = 0, ncol = (motif_size+200), nrow = num_sites)


## Get the signals
cat("Calculating signals...", "\n")
for (a in 1:num_sites){
  cat("site ", a, "of ", num_sites, "\n")
  start <- sites[a,2]-101
  end <- sites[a,3]+100
  width <- end-start
  #
  for (b in 1:width){
    
    # set current bp position on xsome
    pos <- (start+b-1)
    
    # 
    if (pos %in% uq == TRUE){
      sites_signals[a,b] <- ins[(which(ins[,1] == pos)),2]
      
    } # end if (pos %in% uq == TRUE)
  } # end for (b in 1:width)
} # end for (a in 1:num_sites)


## Save
cat("Saving signals...", "\n")
save.image(file = out_filepath)
gc()


