

##
load("~/atac/atac1/mdst8/wt01/footprints/parsed/MDST8-WT-01.OVOL1.motif1.info.Rdata")


signals <- parsedSitesInfo[["combinedbfPassPeakSignal"]][["signal"]]
numsites <- length(signals[,1])
numbp <- length(signals[1,])
motifbp <- numbp-200


## Make a matrix to hold the values
scores <- matrix(data=NA, ncol = 9, nrow = numsites)
colnames(scores) <- c("background", "motif", "motif flank", "total", "flank - motif", "background avg", "flank avg", "flank avg / background avg", "motif avg")

## loop and calculate
for (a in 1:numsites){
  # background
  scores[a,1] <- sum(signals[a,1:90]) + sum(signals[1,(110+motifbp):numbp])
  # motif
  scores[a,2] <- sum(signals[a,100:(100+motifbp)])
  # motif flank
  scores[a,3] <- sum(signals[a,91:100]) + sum(signals[a,(100+motifbp):(110+motifbp)])
  # total
  scores[a,4] <- sum(signals[a,])
  # flank - motif
  scores[a,5] <- scores[a,3] - scores[a,2]
  # background average
  scores[a,6] <- scores[a,1]/180
  # flank average
  scores[a,7] <- scores[a,3]/20
  # flank avg / backgroun avg
  scores[a,8] <- scores[a,7] / scores[a,6]
  # motif avg
  scores[a,9] <- scores[a,2] / motifbp
}