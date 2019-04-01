
##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SRAdb", version = "3.8")

install.packages("rlang")

##
library(SRAdb)
##
dest <- "C:\\Users\\jsk33\\Desktop\\"
##
sraDB <- getSRAdbFile(destdir = dest)
##
sqlfile <- "C:\\Users\\jsk33\\Desktop\\SRAmetadb.sqlite"



sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb.sqlite')

sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)

ac <- "SRX5415026"


getFASTQfile(ac, sra_con, destDir = dest, srcType = "ftp")


dbDisconnect(sra_con)