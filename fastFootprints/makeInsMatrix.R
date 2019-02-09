# Load libraries
gc()
cat("Loading libraries...", "\n")

# Snakemake variables
cat("Loading snakemake variables...", "\n")
counts_path <- snakemake@input[[1]]
out_filepath <- snakemake@output[[1]]


# Load the count matrices
cat("Loading count matrices", "\n")
count_matrix = read.table("/home/ubuntu1/atac/ls1034/wt01/counts/LS1034-WT-01.all.chr1.counts.txt", 
                          sep=" ", 
                          col.names=c("chr", "5", "3"), 
                          fill=FALSE, 
                          strip.white=TRUE)

#
uq <- unique(count_matrix[,2])
# the number of bp with at least 1 insertion
num_uq <- length(uq)

# aggregate the insertions
cat("Aggregating counts...", "\n")
ins <- matrix(data = NA, nrow = num_uq, ncol = 2)
## count them up
for (a in 1:num_uq){
  #cat("Processing site ", a, "of ", num_uq, "\n")
  ins[a,1] <- uq[a]
  ins[a,2] <- length(which(count_matrix[,2] == uq[a]))}

save.image(file = out_filepath)
gc()
cat("Finished.", "\n")
