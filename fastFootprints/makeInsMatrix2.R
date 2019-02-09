
# Snakemake variables
cat("Loading snakemake variables...", "\n")
counts_path <- snakemake@input[[1]]
out_filepath <- snakemake@output[[1]]
# sample_name <- snakemake@wildcards[["sample"]]

# Load the count matrices
cat("Loading count matrices...", "\n")
count_matrix = read.table(counts_path, 
                          sep=" ", 
                          col.names=c("chr", "5", "3"), 
                          fill=FALSE, 
                          strip.white=TRUE)

#
cat("Finding unique values...", "\n")
uq <- unique(count_matrix[,2])
# the number of bp with at least 1 insertion
num_uq <- length(uq)

# aggregate the insertions
cat("Aggregating counts...", "\n")
ins <- matrix(data = NA, nrow = num_uq, ncol = 2)
## count them up
for (a in 1:num_uq){
  cat("Processing site ", a, "of ", num_uq, "\n")
  ins[a,1] <- uq[a]
  ins[a,2] <- length(which(count_matrix[,2] == uq[a]))}

save.image(file = out_filepath)
gc()
cat("Finished.", "\n")
