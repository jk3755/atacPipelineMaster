
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
##
sra_dbname <- sqlfile
##
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)

##
ac <- c("SRS2200162")
##
getFASTQfile(
            in_acc = c("SRR8618987"),
            sra_con = sra_con,
            destDir = dest,
            srcType = "ftp",
            ascpCMD = NULL
            )

getSRAfile(
          in_acc = c("SRX5415026"),
          sra_con = sra_con,
          destDir = dest,
          fileType = 'sra',
          srcType = 'ftp',
          makeDirectory = FALSE,
          method = 'curl',
          ascpCMD = NULL
          )


#############




##
dbDisconnect(sra_con)