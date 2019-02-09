

## Snakemake variables
cat("Loading snakemake variables...", "\n")
signal_path <- snakemake@input[[1]]
## Make a vector with filepaths for all 22 chromosomes
filepaths <- c()
chroms <- paste0("chr", c(1:22))
#
for (a in 1:length(chroms)){
  filepaths[a] <- gsub("chr\\d|chr\\d\\d", chroms[a], signal_path)
}

#
out_filepath <- snakemake@output[[1]]






