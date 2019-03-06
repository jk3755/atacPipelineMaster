library(ATACseqQC)

##
bamfile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\snu61\\wt01\\bam\\SNU61-WT-01-REP1.u.bam"
baifile <- "C:\\Users\\jsk33\\Documents\\atac\\mydata\\snu61\\wt01\\bam\\SNU61-WT-01-REP1.u.bai"

test <- fragSizeDist(
  bamFiles = bamfile,
  bamFiles.labels = "snu61",
  index = baifile
)