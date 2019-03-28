---
  title: "atacFParacne"
author: "Jordan S. Kesner"
date: "January 29, 2019"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/S4Vectors/html/Rle-class.html

```{r}
#install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("mygene", version = "3.8")
#BiocManager::install("chromoMap", version = "3.8")
#
library(ATACseqQC)
library(GenomicRanges)
library(aracne.networks)
library(mygene)
library(chromoMap)
```

Load the required data
```{r}

#
load("C:\\Users\\jsk33\\Desktop\\LS1034-WT-01.CBFA2T2.motif1.Rdata")
l <- ls()
l <- l[-8]
rm(list = l)
rm(l)
#
wg_sites <- FPdata[["wg_bindingsites"]]

#
chromosome <- wg_sites@seqnames
sites <- wg_sites@ranges

# convert the Rle factor object to a linear vector
chromosome <- as.vector(chromosome, mode="any")

```

Step 1 - Build a GRanges object based on the predicted ARACNe targets of current gene
```{r}

# get the ARACNe interactome
coad_interactome <- aracne.networks::reguloncoad

# get the mnx1 targets
mnx1_targets <- names(coad_interactome[["9139"]][["tfmode"]])

# Retrieve info for the network targets
target_list <- list()
for (a in 1:length(mnx1_targets)){
  target_list[a] <- mygene::getGene(geneid = mnx1_targets[a], fields = "all")}

# count the number of gene locations we have
loc_count <- 0
for (a in 1:length(mnx1_targets)){
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      loc_count <- (loc_count + length(target_list[[a]][["genomic_pos"]]))
    } else {
      loc_count <- (loc_count + 1)}}}

# Retrieve the genomic coordinates of all network targets
# retrieve all genomic coords here, can later trim non standard ones more easily with GRanges functions
target_locations <- matrix(data = NA, nrow = loc_count, ncol = 4)
colnames(target_locations) <- c("gene", "chr", "start", "end")
idx <- 1
#
for (a in 1:loc_count){
  
  if (is.null(target_list[[a]][["genomic_pos"]])){next} else {
    
    if (is.list(target_list[[a]][["genomic_pos"]][[1]])){
      
      for (b in 1:length(target_list[[a]][["genomic_pos"]])){
        target_locations[idx,1] <- target_list[[a]][["symbol"]]
        target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][[b]][["chr"]])
        target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][[b]][["start"]]
        target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][[b]][["end"]]
        idx <- (idx+1)
      } # end for (b in 1:length(target_list[[a]][["genomic_pos"]]))
      
    } else {
      
      target_locations[idx,1] <- target_list[[a]][["symbol"]]
      target_locations[idx,2] <- paste0("chr", target_list[[a]][["genomic_pos"]][["chr"]])
      target_locations[idx,3] <- target_list[[a]][["genomic_pos"]][["start"]]
      target_locations[idx,4] <- target_list[[a]][["genomic_pos"]][["end"]]
      idx <- (idx+1)}}}

## Convert the genomic locations of the targets into GRanges
gr <- GRanges(
  seqnames = target_locations[,2],
  ranges = IRanges(start = as.numeric(target_locations[,3]), end = as.numeric(target_locations[,4])))
#strand = c(rep("+", times = num_names)),
#seqinfo = Seqinfo(genome="hg38"),
#score = c(rep(1, times = num_names)))

# prune to standard xsomes
gr <- keepStandardChromosomes(gr, pruning.mode="coarse")


## Intersect the targets GRanges with the binding sites list
intersection <- intersect(gr, wg_sites)

```



Code testing / sanity checks
```{r}

canno <- data.frame(
  name = c(1:148),
  chrom = as.vector(gr@seqnames, mode="any"),
  start = gr@ranges@start,
  data = gr@ranges@width)


cmap <- chromoMap(canno, type = "annotation")
cmap <- chromoMap(gr, type ="heatmap-single", HeatColRange = c("blue","white","red"))



gr1 <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges(start = c(1000, 1000), end = c(2000, 2000))
)

gr2 <- GRanges(
  seqnames = c("chr1", "chr3"),
  ranges = IRanges(start = c(1900, 1000), end = c(3000, 2000))
)

gr3 <- intersect(gr1, gr2)


```




---
  title: "subsetARACNe"
author: "Jordan S. Kesner"
date: "November 5, 2018"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  This script is used to modify a specific ARACNe interactome/regulon for use with VIPER analysis
Gene IDs are in ENTREZ format

Install and load packages
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#BiocManager::install("viper")
#install.packages("rlist")
library(aracne.networks)
library(org.Hs.eg.db)
library(annotate)
library(rlist)
library(viper)
```

Write COAD interactome to text file - this will take a long time  (~30-60 mins)
```{r}
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/") # Write the full interactome to a text file
```

The list of good vs. poor prognosis MRs with motifs:
  TCF7 (6932)
MNX1 (3110)
POU5F1B (5462)
ESRRA (2101)
CDX2 (1045)
HNF4A (3172)
GMEB2 (26205)
HOXA3 (3200)
OVOL1 (5017)
ASCL2 (430)
ZSWIM1 (90204)
CBFA2T2 (9139)

FIRE
ADNP (23394)
ZNF696 (79943)
ZMYND8 (23613)
TAF4 (6874)

Write the MR regulons to file
```{r}
data(reguloncoad)
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/TCF7.txt", regulator="6932")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/MNX1.txt", regulator="3110")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/POU5F1B.txt", regulator="5462")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ESRRA.txt", regulator="2101")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/CDX2.txt", regulator="1045")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/HNF4A.txt", regulator="3172")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/GMEB2.txt", regulator="26205")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/HOXA3.txt", regulator="3200")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/OVOL1.txt", regulator="5017")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ASCL2.txt", regulator="430")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZSWIM1.txt", regulator="90204")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/CBFA2T2.txt", regulator="9139")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ADNP.txt", regulator="23394")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZNF696.txt", regulator="79943")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZMYND8.txt", regulator="23613")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/TAF4.txt", regulator="6874")
```

Read in the saved regulons
```{r}
## Note, must delete the first line in the text files first

# TCF7
coad_tcf7 <- matrix(data=NA, nrow=71, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/TCF7.txt"
con <- file(description=file,open="r")
for(i in 1:70){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_tcf7[i,1] <- tmp[1]
  coad_tcf7[i,2] <- tmp[2]
  coad_tcf7[i,3] <- tmp[3]
  coad_tcf7[i,4] <- tmp[4]}
close(con)


# MNX1
coad_mnx1 <- matrix(data=NA, nrow=107, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/MNX1.txt"
con <- file(description=file,open="r")
for(i in 1:107){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_mnx1[i,1] <- tmp[1]
  coad_mnx1[i,2] <- tmp[2]
  coad_mnx1[i,3] <- tmp[3]
  coad_mnx1[i,4] <- tmp[4]}
close(con)

# POU5F1B
coad_pou5f1b <- matrix(data=NA, nrow=227, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/POU5F1B.txt"
con <- file(description=file,open="r")
for(i in 1:227){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_pou5f1b[i,1] <- tmp[1]
  coad_pou5f1b[i,2] <- tmp[2]
  coad_pou5f1b[i,3] <- tmp[3]
  coad_pou5f1b[i,4] <- tmp[4]}
close(con)

# ESRRA
coad_esrra <- matrix(data=NA, nrow=151, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ESRRA.txt"
con <- file(description=file,open="r")
for(i in 1:151){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_esrra[i,1] <- tmp[1]
  coad_esrra[i,2] <- tmp[2]
  coad_esrra[i,3] <- tmp[3]
  coad_esrra[i,4] <- tmp[4]}
close(con)

# CDX2
coad_cdx2 <- matrix(data=NA, nrow=97, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/CDX2.txt"
con <- file(description=file,open="r")
for(i in 1:97){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_cdx2[i,1] <- tmp[1]
  coad_cdx2[i,2] <- tmp[2]
  coad_cdx2[i,3] <- tmp[3]
  coad_cdx2[i,4] <- tmp[4]}
close(con)

# HNF4A
coad_hnf4a <- matrix(data=NA, nrow=91, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/HNF4A.txt"
con <- file(description=file,open="r")
for(i in 1:91){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_hnf4a[i,1] <- tmp[1]
  coad_hnf4a[i,2] <- tmp[2]
  coad_hnf4a[i,3] <- tmp[3]
  coad_hnf4a[i,4] <- tmp[4]}
close(con)

# GMEB2
coad_gmeb2 <- matrix(data=NA, nrow=289, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/GMEB2.txt"
con <- file(description=file,open="r")
for(i in 1:289){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_gmeb2[i,1] <- tmp[1]
  coad_gmeb2[i,2] <- tmp[2]
  coad_gmeb2[i,3] <- tmp[3]
  coad_gmeb2[i,4] <- tmp[4]}
close(con)

# HOXA3
coad_hoxa3 <- matrix(data=NA, nrow=93, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/HOXA3.txt"
con <- file(description=file,open="r")
for(i in 1:93){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_hoxa3[i,1] <- tmp[1]
  coad_hoxa3[i,2] <- tmp[2]
  coad_hoxa3[i,3] <- tmp[3]
  coad_hoxa3[i,4] <- tmp[4]}
close(con)

# OVOL1
coad_ovol1 <- matrix(data=NA, nrow=149, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/OVOL1.txt"
con <- file(description=file,open="r")
for(i in 1:149){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_ovol1[i,1] <- tmp[1]
  coad_ovol1[i,2] <- tmp[2]
  coad_ovol1[i,3] <- tmp[3]
  coad_ovol1[i,4] <- tmp[4]}
close(con)

# ASCL2
coad_ascl2 <- matrix(data=NA, nrow=127, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ASCL2.txt"
con <- file(description=file,open="r")
for(i in 1:127){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_ascl2[i,1] <- tmp[1]
  coad_ascl2[i,2] <- tmp[2]
  coad_ascl2[i,3] <- tmp[3]
  coad_ascl2[i,4] <- tmp[4]}
close(con)

# ZSWIM1
coad_zswim1 <- matrix(data=NA, nrow=108, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZSWIM1.txt"
con <- file(description=file,open="r")
for(i in 1:108){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_zswim1[i,1] <- tmp[1]
  coad_zswim1[i,2] <- tmp[2]
  coad_zswim1[i,3] <- tmp[3]
  coad_zswim1[i,4] <- tmp[4]}
close(con)

# CBFA2T2
coad_cbfa2t2 <- matrix(data=NA, nrow=152, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/CBFA2T2.txt"
con <- file(description=file,open="r")
for(i in 1:152){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_cbfa2t2[i,1] <- tmp[1]
  coad_cbfa2t2[i,2] <- tmp[2]
  coad_cbfa2t2[i,3] <- tmp[3]
  coad_cbfa2t2[i,4] <- tmp[4]}

# ADNP
coad_adnp <- matrix(data=NA, nrow=146, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ADNP.txt"
con <- file(description=file,open="r")
for(i in 1:146){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_adnp[i,1] <- tmp[1]
  coad_adnp[i,2] <- tmp[2]
  coad_adnp[i,3] <- tmp[3]
  coad_adnp[i,4] <- tmp[4]}
close(con)

# ZNF696
coad_znf696 <- matrix(data=NA, nrow=74, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZNF696.txt"
con <- file(description=file,open="r")
for(i in 1:74){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_znf696[i,1] <- tmp[1]
  coad_znf696[i,2] <- tmp[2]
  coad_znf696[i,3] <- tmp[3]
  coad_znf696[i,4] <- tmp[4]}

# ZMYND8
coad_zmynd8 <- matrix(data=NA, nrow=127, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZMYND8.txt"
con <- file(description=file,open="r")
for(i in 1:127){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_zmynd8[i,1] <- tmp[1]
  coad_zmynd8[i,2] <- tmp[2]
  coad_zmynd8[i,3] <- tmp[3]
  coad_zmynd8[i,4] <- tmp[4]}
close(con)

# TAF4
coad_taf4 <- matrix(data=NA, nrow=122, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/TAF4.txt"
con <- file(description=file,open="r")
for(i in 1:122){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_taf4[i,1] <- tmp[1]
  coad_taf4[i,2] <- tmp[2]
  coad_taf4[i,3] <- tmp[3]
  coad_taf4[i,4] <- tmp[4]}
close(con)

```

Make a list of the regulons
```{r}
coad_regulons <- list()
coad_regulons$adnp <- coad_adnp
coad_regulons$ascl2 <- coad_ascl2
coad_regulons$cbfa2t2 <- coad_cbfa2t2
coad_regulons$cdx2 <- coad_cdx2
coad_regulons$esrra <- coad_esrra
coad_regulons$gmeb2 <- coad_gmeb2
coad_regulons$hnf4a <- coad_hnf4a
coad_regulons$hoxa3 <- coad_hoxa3
coad_regulons$mnx1 <- coad_mnx1
coad_regulons$ovol1 <- coad_ovol1
coad_regulons$pou5f1b <- coad_pou5f1b
coad_regulons$taf4 <- coad_taf4
coad_regulons$tcf7 <- coad_tcf7
coad_regulons$zmynd8 <- coad_zmynd8
coad_regulons$znf696 <- coad_znf696
coad_regulons$zswim1 <- coad_zswim1

# remove the single versions, if wanted
rm(coad_adnp, coad_ascl2, coad_cbfa2t2, coad_cdx2, coad_esrra, coad_gmeb2, coad_hnf4a, coad_hoxa3, coad_mnx1, coad_ovol1, coad_pou5f1b, coad_taf4, coad_tcf7, coad_zmynd8, coad_znf696, coad_zswim1)

```





Convert between ENTREZ and GENE SYMBOL ** working on it **
  ```{r}

```









Code to tinker with interactome
```{r}
gene_id <- "10002"
setwd("C:/Users/Jordan/Desktop/")
#interactome_tinker <- read.table(file = "COAD_interactome.txt", sep = "\t") # read in the interactome from text file
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1] == gene_id),] # getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1] == gene_id)),] # removing the original regulon from interactome
#colnames(regulon_motifMI) <- colnames(regulon_req)
write.table(interactome_tinker, file = "COAD_interactome_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE) # write the modified interactome to file

mod_interactome <- aracne2regulon("COAD_interactome_tinkered.txt", dset, format = "3col", verbose = TRUE) # regenerate the new interactome

print(regul)
#coadViper <- viper(dset, regul, method = "scale") #the viper algorithm
#dim(coadViper)


### below is directly from Ajay
#2018 Nov 6
#code to tinker with interactome
library(viper)
setwd("C:/Users/Jordan/Desktop/")
# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")
load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
geneId <- "4089" #smad4
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1]==geneId),]#getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1]==geneId)),]# removing the req regulon from interactome
colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
write.table(interactome_tinker, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)
write.table(regulon_req, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
regul <- aracne2regulon("network_coad_tfCotf_3col_tinkered.txt", dset, format = "3col", verbose = TRUE)#regulon
#
#
print(regul)

coadViper <- viper(dset, regul, method = "scale") # the viper algorithm. takes only a few minutes to run
dim(coadViper)

plot(coadViper, mrs = 10, color = c("cornflowerblue", "salmon"),
     pval = NULL, bins = 500, cex = 0, density = 0, smooth = 0,
     sep = 0.2, hybrid = TRUE, include = c("expression", "activity"),
     gama = 2)



```

removeTargets
```{r}
## A function to remove targets from all regulator regulons in network
## This is used when a specific set of genes is determined to be in a region of constituitive heterochromatin
## Not sure if it will work, but can try to set tfmode to 0 instead of removing edge
removeTargets <- function(interactome, genes){
  
}

a <- c(as.character("158056"), as.character("57096"))
list.remove(reguloncoad[["10002"]], a)
lis <- vector("list", 3)
lis[[3]] <- NULL

list.remove(reguloncoad[["10009"]][["tfmode"]][["10735"]], c("10735"))


names <- as.character(names(reguloncoad)) # get the entrez ids in the list
symbols <- c()
for (a in 1:length(names)){symbols[a] <- getSYMBOL(names[a], "org.Hs.eg")} # convert the entrez ids to gene symbols


## Remove a set of target genes from all regulator regulons
## This should be used when a gene is determined to be located in constituitive heterochromatin
remove_list <- c("TCF7")
remove_num <- length(remove_list)

reguloncoad[["10002"]]

# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")
#load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
#geneId <- "4089" #smad4

data(package="aracne.networks")$results[, "Item"] # show available packages
data(regulonblca) 
write.regulon(regulonblca, n = 10) # write regulon to file, with 10 interactions
data(regulonblca)
write.regulon(regulonblca, regulator="399") # analyze a specific regulon
getSYMBOL("6932", "org.Hs.eg") # this function converts an ENTREZ ID into a gene symbol

data(reguloncoad) # load the interactome from the package
reguloncoad # data is loaded as a promise object, must interact with it to conver to df
regulators <- names(reguloncoad) # returns list of all regulators in interactome
targets <- c() # a vector to hold the targets
for (a in 1:length(regulators)){targets <- c(targets, names(reguloncoad[[a]][["tfmode"]]))} # populate targets vector
unique <- c(regulators, targets) # merge targets and regulators list
unique <- unique(unique) # get unique entries
gene_id <- matrix(data = NA, ncol = 2, nrow = length(unique)) # a matrix that will hold ENTREZ and SYMBOL for all unique genes
for (b in 1:length(unique)){
  ent <- as.character(unique[[b]]); sym <- as.character(getSYMBOL(ent, "org.Hs.eg"))
  gene_id[b,1] <- ent; gene_id[b,2] <- sym}
rm(regulators, targets, a, b, ent, sym, unique) # cleanup

colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
```

---
  title: "subsetARACNe"
author: "Jordan S. Kesner"
date: "November 5, 2018"
output: html_document
editor_options: 
  chunk_output_type: console
---
  Code for working with ARACNe and ATACseq

This script is used to modify a specific ARACNe interactome/regulon for use with VIPER analysis
Gene IDs are in ENTREZ format

Install and load packages
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
#BiocManager::install("viper")
#install.packages("rlist")
library(aracne.networks)
library(org.Hs.eg.db)
library(annotate)
library(rlist)
library(viper)
```

When given a specific ARACNe interactome, this script will scan the object for all regulators, translate them to gene symbols, and format a string to be plugged in directly to the ATACseqFP pipeline as a gene list. This will allow you to do a full scan of all predicted regulators for a given cell type for TF footprints

Run operation
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("aracne.networks", version = "3.8")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("annotate")
library(aracne.networks)
library(org.Hs.eg.db)
library(annotate)

regulators <- names(reguloncoad) # returns list of all regulators in interactome
symbols <- c()
for (a in 1:length(regulators)){
  symbols[a] <- getSYMBOL(regulators[a],"org.Hs.eg")
}
# The symbols vector can be piped directly into the pipeline now
```






Write COAD interactome to text file - this will take a long time  (~30-60 mins)
```{r}
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/") # Write the full interactome to a text file
```

The list of good vs. poor prognosis MRs with motifs:
  TCF7 (6932)
MNX1 (3110)
POU5F1B (5462)
ESRRA (2101)
CDX2 (1045)
HNF4A (3172)
GMEB2 (26205)
HOXA3 (3200)
OVOL1 (5017)
ASCL2 (430)
ZSWIM1 (90204)
CBFA2T2 (9139)

FIRE
ADNP (23394)
ZNF696 (79943)
ZMYND8 (23613)
TAF4 (6874)

Write the MR regulons to file
```{r}
data(reguloncoad)
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/TCF7.txt", regulator="6932")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/MNX1.txt", regulator="3110")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/POU5F1B.txt", regulator="5462")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ESRRA.txt", regulator="2101")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/CDX2.txt", regulator="1045")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/HNF4A.txt", regulator="3172")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/GMEB2.txt", regulator="26205")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/HOXA3.txt", regulator="3200")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/OVOL1.txt", regulator="5017")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ASCL2.txt", regulator="430")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZSWIM1.txt", regulator="90204")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/CBFA2T2.txt", regulator="9139")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ADNP.txt", regulator="23394")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZNF696.txt", regulator="79943")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/ZMYND8.txt", regulator="23613")
write.regulon(reguloncoad, file = "C:/Users/Jordan/Desktop/COAD_regulons/TAF4.txt", regulator="6874")
```

Read in the saved regulons
```{r}
## Note, must delete the first line in the text files first

# TCF7
coad_tcf7 <- matrix(data=NA, nrow=71, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/TCF7.txt"
con <- file(description=file,open="r")
for(i in 1:70){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_tcf7[i,1] <- tmp[1]
  coad_tcf7[i,2] <- tmp[2]
  coad_tcf7[i,3] <- tmp[3]
  coad_tcf7[i,4] <- tmp[4]}
close(con)


# MNX1
coad_mnx1 <- matrix(data=NA, nrow=107, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/MNX1.txt"
con <- file(description=file,open="r")
for(i in 1:107){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_mnx1[i,1] <- tmp[1]
  coad_mnx1[i,2] <- tmp[2]
  coad_mnx1[i,3] <- tmp[3]
  coad_mnx1[i,4] <- tmp[4]}
close(con)

# POU5F1B
coad_pou5f1b <- matrix(data=NA, nrow=227, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/POU5F1B.txt"
con <- file(description=file,open="r")
for(i in 1:227){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_pou5f1b[i,1] <- tmp[1]
  coad_pou5f1b[i,2] <- tmp[2]
  coad_pou5f1b[i,3] <- tmp[3]
  coad_pou5f1b[i,4] <- tmp[4]}
close(con)

# ESRRA
coad_esrra <- matrix(data=NA, nrow=151, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ESRRA.txt"
con <- file(description=file,open="r")
for(i in 1:151){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_esrra[i,1] <- tmp[1]
  coad_esrra[i,2] <- tmp[2]
  coad_esrra[i,3] <- tmp[3]
  coad_esrra[i,4] <- tmp[4]}
close(con)

# CDX2
coad_cdx2 <- matrix(data=NA, nrow=97, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/CDX2.txt"
con <- file(description=file,open="r")
for(i in 1:97){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_cdx2[i,1] <- tmp[1]
  coad_cdx2[i,2] <- tmp[2]
  coad_cdx2[i,3] <- tmp[3]
  coad_cdx2[i,4] <- tmp[4]}
close(con)

# HNF4A
coad_hnf4a <- matrix(data=NA, nrow=91, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/HNF4A.txt"
con <- file(description=file,open="r")
for(i in 1:91){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_hnf4a[i,1] <- tmp[1]
  coad_hnf4a[i,2] <- tmp[2]
  coad_hnf4a[i,3] <- tmp[3]
  coad_hnf4a[i,4] <- tmp[4]}
close(con)

# GMEB2
coad_gmeb2 <- matrix(data=NA, nrow=289, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/GMEB2.txt"
con <- file(description=file,open="r")
for(i in 1:289){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_gmeb2[i,1] <- tmp[1]
  coad_gmeb2[i,2] <- tmp[2]
  coad_gmeb2[i,3] <- tmp[3]
  coad_gmeb2[i,4] <- tmp[4]}
close(con)

# HOXA3
coad_hoxa3 <- matrix(data=NA, nrow=93, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/HOXA3.txt"
con <- file(description=file,open="r")
for(i in 1:93){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_hoxa3[i,1] <- tmp[1]
  coad_hoxa3[i,2] <- tmp[2]
  coad_hoxa3[i,3] <- tmp[3]
  coad_hoxa3[i,4] <- tmp[4]}
close(con)

# OVOL1
coad_ovol1 <- matrix(data=NA, nrow=149, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/OVOL1.txt"
con <- file(description=file,open="r")
for(i in 1:149){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_ovol1[i,1] <- tmp[1]
  coad_ovol1[i,2] <- tmp[2]
  coad_ovol1[i,3] <- tmp[3]
  coad_ovol1[i,4] <- tmp[4]}
close(con)

# ASCL2
coad_ascl2 <- matrix(data=NA, nrow=127, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ASCL2.txt"
con <- file(description=file,open="r")
for(i in 1:127){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_ascl2[i,1] <- tmp[1]
  coad_ascl2[i,2] <- tmp[2]
  coad_ascl2[i,3] <- tmp[3]
  coad_ascl2[i,4] <- tmp[4]}
close(con)

# ZSWIM1
coad_zswim1 <- matrix(data=NA, nrow=108, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZSWIM1.txt"
con <- file(description=file,open="r")
for(i in 1:108){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_zswim1[i,1] <- tmp[1]
  coad_zswim1[i,2] <- tmp[2]
  coad_zswim1[i,3] <- tmp[3]
  coad_zswim1[i,4] <- tmp[4]}
close(con)

# CBFA2T2
coad_cbfa2t2 <- matrix(data=NA, nrow=152, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/CBFA2T2.txt"
con <- file(description=file,open="r")
for(i in 1:152){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_cbfa2t2[i,1] <- tmp[1]
  coad_cbfa2t2[i,2] <- tmp[2]
  coad_cbfa2t2[i,3] <- tmp[3]
  coad_cbfa2t2[i,4] <- tmp[4]}

# ADNP
coad_adnp <- matrix(data=NA, nrow=146, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ADNP.txt"
con <- file(description=file,open="r")
for(i in 1:146){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_adnp[i,1] <- tmp[1]
  coad_adnp[i,2] <- tmp[2]
  coad_adnp[i,3] <- tmp[3]
  coad_adnp[i,4] <- tmp[4]}
close(con)

# ZNF696
coad_znf696 <- matrix(data=NA, nrow=74, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZNF696.txt"
con <- file(description=file,open="r")
for(i in 1:74){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_znf696[i,1] <- tmp[1]
  coad_znf696[i,2] <- tmp[2]
  coad_znf696[i,3] <- tmp[3]
  coad_znf696[i,4] <- tmp[4]}

# ZMYND8
coad_zmynd8 <- matrix(data=NA, nrow=127, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/ZMYND8.txt"
con <- file(description=file,open="r")
for(i in 1:127){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_zmynd8[i,1] <- tmp[1]
  coad_zmynd8[i,2] <- tmp[2]
  coad_zmynd8[i,3] <- tmp[3]
  coad_zmynd8[i,4] <- tmp[4]}
close(con)

# TAF4
coad_taf4 <- matrix(data=NA, nrow=122, ncol=4)
file <- "C:/Users/Jordan/Desktop/COAD_regulons/TAF4.txt"
con <- file(description=file,open="r")
for(i in 1:122){
  tmp <- scan(file=con, nlines = 1, quiet=TRUE)
  coad_taf4[i,1] <- tmp[1]
  coad_taf4[i,2] <- tmp[2]
  coad_taf4[i,3] <- tmp[3]
  coad_taf4[i,4] <- tmp[4]}
close(con)

```

Make a list of the regulons
```{r}
coad_regulons <- list()
coad_regulons$adnp <- coad_adnp
coad_regulons$ascl2 <- coad_ascl2
coad_regulons$cbfa2t2 <- coad_cbfa2t2
coad_regulons$cdx2 <- coad_cdx2
coad_regulons$esrra <- coad_esrra
coad_regulons$gmeb2 <- coad_gmeb2
coad_regulons$hnf4a <- coad_hnf4a
coad_regulons$hoxa3 <- coad_hoxa3
coad_regulons$mnx1 <- coad_mnx1
coad_regulons$ovol1 <- coad_ovol1
coad_regulons$pou5f1b <- coad_pou5f1b
coad_regulons$taf4 <- coad_taf4
coad_regulons$tcf7 <- coad_tcf7
coad_regulons$zmynd8 <- coad_zmynd8
coad_regulons$znf696 <- coad_znf696
coad_regulons$zswim1 <- coad_zswim1

# remove the single versions, if wanted
rm(coad_adnp, coad_ascl2, coad_cbfa2t2, coad_cdx2, coad_esrra, coad_gmeb2, coad_hnf4a, coad_hoxa3, coad_mnx1, coad_ovol1, coad_pou5f1b, coad_taf4, coad_tcf7, coad_zmynd8, coad_znf696, coad_zswim1)

```





Convert between ENTREZ and GENE SYMBOL ** working on it **
  ```{r}

```









Code to tinker with interactome
```{r}
gene_id <- "10002"
setwd("C:/Users/Jordan/Desktop/")
#interactome_tinker <- read.table(file = "COAD_interactome.txt", sep = "\t") # read in the interactome from text file
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1] == gene_id),] # getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1] == gene_id)),] # removing the original regulon from interactome
#colnames(regulon_motifMI) <- colnames(regulon_req)
write.table(interactome_tinker, file = "COAD_interactome_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE) # write the modified interactome to file

mod_interactome <- aracne2regulon("COAD_interactome_tinkered.txt", dset, format = "3col", verbose = TRUE) # regenerate the new interactome

print(regul)
#coadViper <- viper(dset, regul, method = "scale") #the viper algorithm
#dim(coadViper)


### below is directly from Ajay
#2018 Nov 6
#code to tinker with interactome
library(viper)
setwd("C:/Users/Jordan/Desktop/")
# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")
load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
geneId <- "4089" #smad4
interactome_tinker <- read.table(file = "network_coad_tfCotf_3col_sorted.txt", sep = "\t")
regulon_req <- interactome_tinker[which(interactome_tinker[,1]==geneId),]#getting req regulon to tinker
interactome_tinker <- interactome_tinker[-(which(interactome_tinker[,1]==geneId)),]# removing the req regulon from interactome
colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
write.table(interactome_tinker, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = FALSE, row.names = FALSE, col.names = FALSE)
write.table(regulon_req, file = "network_coad_tfCotf_3col_tinkered.txt", sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
regul <- aracne2regulon("network_coad_tfCotf_3col_tinkered.txt", dset, format = "3col", verbose = TRUE)#regulon
#
#
print(regul)

coadViper <- viper(dset, regul, method = "scale") # the viper algorithm. takes only a few minutes to run
dim(coadViper)

plot(coadViper, mrs = 10, color = c("cornflowerblue", "salmon"),
     pval = NULL, bins = 500, cex = 0, density = 0, smooth = 0,
     sep = 0.2, hybrid = TRUE, include = c("expression", "activity"),
     gama = 2)



```

removeTargets
```{r}
## A function to remove targets from all regulator regulons in network
## This is used when a specific set of genes is determined to be in a region of constituitive heterochromatin
## Not sure if it will work, but can try to set tfmode to 0 instead of removing edge
removeTargets <- function(interactome, genes){
  
}

a <- c(as.character("158056"), as.character("57096"))
list.remove(reguloncoad[["10002"]], a)
lis <- vector("list", 3)
lis[[3]] <- NULL

list.remove(reguloncoad[["10009"]][["tfmode"]][["10735"]], c("10735"))


names <- as.character(names(reguloncoad)) # get the entrez ids in the list
symbols <- c()
for (a in 1:length(names)){symbols[a] <- getSYMBOL(names[a], "org.Hs.eg")} # convert the entrez ids to gene symbols


## Remove a set of target genes from all regulator regulons
## This should be used when a gene is determined to be located in constituitive heterochromatin
remove_list <- c("TCF7")
remove_num <- length(remove_list)

reguloncoad[["10002"]]

# save(dset,regulon_motifMI, file = "dset_motifMIregulon.RData")
#load(file = "dset_motifMIregulon.RData")#get the dset and regulon_motifMI
#geneId <- "4089" #smad4

data(package="aracne.networks")$results[, "Item"] # show available packages
data(regulonblca) 
write.regulon(regulonblca, n = 10) # write regulon to file, with 10 interactions
data(regulonblca)
write.regulon(regulonblca, regulator="399") # analyze a specific regulon
getSYMBOL("6932", "org.Hs.eg") # this function converts an ENTREZ ID into a gene symbol

data(reguloncoad) # load the interactome from the package
reguloncoad # data is loaded as a promise object, must interact with it to conver to df
regulators <- names(reguloncoad) # returns list of all regulators in interactome
targets <- c() # a vector to hold the targets
for (a in 1:length(regulators)){targets <- c(targets, names(reguloncoad[[a]][["tfmode"]]))} # populate targets vector
unique <- c(regulators, targets) # merge targets and regulators list
unique <- unique(unique) # get unique entries
gene_id <- matrix(data = NA, ncol = 2, nrow = length(unique)) # a matrix that will hold ENTREZ and SYMBOL for all unique genes
for (b in 1:length(unique)){
  ent <- as.character(unique[[b]]); sym <- as.character(getSYMBOL(ent, "org.Hs.eg"))
  gene_id[b,1] <- ent; gene_id[b,2] <- sym}
rm(regulators, targets, a, b, ent, sym, unique) # cleanup

colnames(regulon_motifMI) <- colnames(regulon_req)
regulon_req <- rbind(regulon_req,regulon_motifMI)
regulon_req <- regulon_req[order(regulon_req[,3], decreasing = TRUE),]
```

