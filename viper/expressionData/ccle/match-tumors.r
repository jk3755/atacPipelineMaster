# Match cell-lines to tcga tumors

# Initializing
library(viper)
setwd("~/scratch/darwin/database/ccle/")

# Loading TCGA data
load("../tcga/tcga-vipermat.rda")
vipermat.tcga <- vipermat
load("../tcga/tcga-samples.rda")
samples.tcga <- samples

# Loading ccle data
load("ccle-vipermat.rda")
vipermat.ccle <- vipermat
samples.ccle <- clines

# Constructing the regulons with the top 50 MRs
reg <- apply(vipermat.tcga[rownames(vipermat.tcga) %in% rownames(vipermat.ccle), ], 2, function(x) {
    tfmode <- sign(x[order(x)[c(1:25, (length(x)-24):length(x))]])
    list(tfmode=tfmode, likelihood=rep(1, length(tfmode)))
})
class(reg) <- "regulon"

# Computing the similarity
nes <- aREA(vipermat.ccle, reg)$nes

# checking for tissue type
clines <- clines[clines[, 3]!="-", ]
nes <- nes[, match(clines[, 1], colnames(nes))]

cline.tcga <- lapply(1:ncol(nes), function(i, nes, tt, nn) {
    tt[order(nes[, i], decreasing=TRUE)[1:nn]]
}, tt=samples.tcga[, 1], nes=nes, nn=10)

tmp <- sapply(1:length(cline.tcga), function(i, cline.tcga, tt) {
    tt[i] %in% cline.tcga[[i]]
}, tt=clines[, 3], cline.tcga=cline.tcga)
