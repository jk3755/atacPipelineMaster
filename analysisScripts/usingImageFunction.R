
# image(x, y, z, zlim, xlim, ylim, col = heat.colors(12),
# add = FALSE, xaxs = "i", yaxs = "i", xlab, ylab,
# breaks, oldstyle = FALSE, useRaster, .)

# x, y - locations of grid lines at which the values in z are measured
# z - a numeric or logical matrix containing the values to be plotted
# zlim - minimum and maximum z values for which color should be plotted
# xlim, ylim - limits to the x and y axis
# col - a list of colors (rainbow, heat.colors, topo.colors, terrain.colors)
# xlab, ylab - labels of the x and y axis
# breaks - a set of finite numeric breakpoints for the colors
# other graphical parameters can be passed with ...


# Basic usage
image(sig[1:100,])

# rows are plotted on the x-axis
# columns are plotted on the y-axis
image(
      sig[1:10,],
      col = rainbow(5),
      xlab = "X axis",
      ylab = "Y axis"
      )

## Using Rle lists and ChipPeakAnno, etc
plus <- sigs[["+"]]
minus <- sigs[["-"]]

vplus <- c(plus)
vminus <- c(minus)

l <- NumericList(vplus)

rplus <- RleList(vplus)
rplus@partitioning@end <- rep(207, 1667)
rplus@partitioning@NAMES <- "ko"

## get coverage
cvglist <- sapply(bamIn, coverage)
cvglist <- cvglist[c("+", "-")]
cvglist <- lapply(cvglist, function(.ele)
  .ele[names(.ele) %in% seqlev])


ChIPpeakAnno::featureAlignedSignal(
                                    cvglists = rplus,
                                    feature.gr = sites,
                                    upstream = 100,
                                    downstream = 100,
                                    n.tile = 100
                                  )


# cvglists - Output of featureAlignedSignal or a list of SimpleRleList or RleList
# feature.gr - An object of GRanges with identical width. If the width equal to 1, you can use upstream and downstream to set the range for plot. If the width not equal to 1, you can use zeroAt to set the zero point of the heatmap.
# upstream, downstream - upstream or dwonstream from the feature.gr. It must keep same as featureAlignedSignal. It is used for x-axis label.
# zeroAt - zero point position of feature.gr
# n.tile - The number of tiles to generate for each element of feature.gr, default is 100
# annoMcols - The columns of metadata of feature.gr that specifies the annotations shown of the right side of the heatmap
# sortBy - Sort the feature.gr by columns by annoMcols and then the signals of the given samples. Default is the first sample. Set to NULL to disable sort.
# colors - vector of colors used in heatmap
# lower.extreme, upper.extreme - The lower and upper boundary value of each samples
# margin - margin for the plot region
# gap - gap in between each heatmap column
# newpage - call grid.newpage or not, default = true
# gp - a gpar object that can be used for text

#svg(file = savepath_heatmap)

## USE THIS STRUCTURE FOR HEATMAPS
## first, combine the signals from plus and minus strand
sigs <- parsedSitesInfo[["bfPassPeakSignals"]]
sig_plus <- sigs["+"]
sig_minus <- sigs["-"]
#
numsites <- length(sig_plus[["+"]][,1])
numbp <- length(sig_plus[["+"]][1,])
#
combined <- list()
combined$signal <- matrix(data = NA, nrow = numsites, ncol = numbp)
# combine the signals
for (i in 1:numsites){combined$signal[i,] <- sigs[["+"]][i,] + sigs[["-"]][i,]}

sites <- parsedSitesInfo[["bfPassPeakSites"]]
num_bp <- length(sigs[["+"]][1,]) 
ChIPpeakAnno::featureAlignedHeatmap(combined, 
                                    feature.gr=reCenterPeaks(sites,width=num_bp), 
                                    annoMcols="score",
                                    sortBy="score",
                                    n.tile=num_bp)
#dev.off()


ChIPpeakAnno::featureAlignedHeatmap(
                                    cvglists = sigs$`+`,
                                    feature.gr = sites
                                  )






