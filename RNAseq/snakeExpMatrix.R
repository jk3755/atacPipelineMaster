
## Install/load
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport", version = "3.8")
##
library(tximport)
library(ensembldb)

###
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
Tx <- transcripts(edb,
                  columns = c(listColumns(edb , "tx"), "gene_name"),
                  return.type = "DataFrame")
###
tx2gene <- data.frame(TXNAME = Tx@listData[["tx_id"]], GENEID = Tx@listData[["gene_name"]])
### Import from salmon data output
files <- "C:\\Users\\jsk33\\Desktop\\SRR8618987.snu61.salmon.quant\\quant.sf"
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
head(txi$counts)

###
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)


###
txi.sum <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion = TRUE)
all.equal(txi$counts, txi.sum$counts)
