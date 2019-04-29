## Make a plot for the bf passing sites
#svgPath <- paste0(dirPath, "footprints/graphs/bf/", sampleName, ".", geneName, ".", "motif", a, ".bf.sites.svg")
#svg(file = svgPath)
#cat("Saving peaks footprint image at path:", svgPath, "\n")
#plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".bfsites")
#plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = bfVector)
#dev.off()

# ## HEATMAPS ## CURRENTLY OFF ##
# #### Make heatmap for bf passing sites ####
# ## USE THIS STRUCTURE FOR HEATMAPS
# ## first, combine the signals from plus and minus strand
# heatSigs <- bfInsMatrix
# heatNumSites <- bfNumSites
# heatNumBP <- siteWidth
# heatSites <- bfSites
# ## scale each row individually
# for (f in 1:heatNumSites){
#   maxsig <- max(heatSigs[f,])
#   for (g in 1:heatNumBP){heatSigs[f,g] <- (heatSigs[f,g] / maxsig)}}
# maxsig <- 1
# ## invert signals
# for (h in 1:heatNumSites){for (i in 1:heatNumBP){heatSigs[h,i] <- (1-heatSigs[h,i])}}
# 
# ## Annotate the combined sublist name which will become the tital of the heatmap plot
# heatTitle <- paste0(geneName, "_motif", a, "_numsites", heatNumSites)
# combined <- list()
# com <- paste0("combined$", heatTitle, " <- heatSigs")
# eval(parse(text = com))
# svgPath <- paste0(dirPath, "footprints/graphs/heatmaps/", sampleName, ".", geneName, ".", "motif", a, ".bfpeak.sites.heatmap.svg")
# svg(file = svgPath)
# cat("Saving svg footprint image at path:", svgPath, "\n")
# ## Margin controls
# # margin(a,b,c,d)
# # a = size of graph from top to bottom, higher value = smaller. default = 0.1
# # b = size of graph from left to right, higher value = smaller. default = 0.005
# # c = flips x axis?
# # d = margin from right side of page, higher = smaller. set at 0.2 so legends dont overlap
# # good settings for ATACseq = c(0.1,0.005,0.05,0.2)
# # bias setting >1 puts more colors at higher values, very useful for dealing with washout of low values
# ChIPpeakAnno::featureAlignedHeatmap(combined,
#                                     feature.gr = reCenterPeaks(heatSites, width = heatNumBP),
#                                     upper.extreme = maxsig, # set this to control the heatmap scale
#                                     annoMcols = "score",
#                                     sortBy = "score",
#                                     n.tile = heatNumBP,
#                                     margin = c(0.1, 0.005, 0.05, 0.2),
#                                     color = colorRampPalette(c("white","grey98","grey97","grey99", "firebrick"), bias = 0.9)(100),
#                                     gp = gpar(fontsize = 10),
#                                     newpage = TRUE)
# dev.off()



plotInsProb <- function(plotTitle = c(""), motifWidth, motifPWM, plotLogo = FALSE, insVector){
  
  ## This function uses code adapted from the R package ATACSeqQC
  ## Plot the figure in a new page in the viewport
  grid.newpage()
  
  ## Data
  totalBP <- length(insVector)
  flankBP <- ((totalBP - motifWidth) / 2 ) ## The number of BP flanking the motif on each side
  
  ## Plot information
  xlab = "Dist. to motif (bp)"
  ylab = "Tn5 fragmentation probability"
  xlim <- c(0, totalBP + 1)
  ylim <- c(0, max(insVector) * 1.12)
  
  ## Add the plotting margins to the viewport (sets the outer bounds of the entire image)
  vp <- plotViewport(margins=c(5.1, 5.1, 4.1, 2.1), name="plotRegion")
  pushViewport(vp)
  
  ## Viewport for the graph plotting area
  vp1 <- viewport(y=.4, height=.8,
                  xscale=xlim,
                  yscale=ylim,
                  name="footprints")
  pushViewport(vp1)
  
  ## Add the insertion probability data line
  grid.lines(x=1:totalBP,
             y=insVector,
             default.units="native",
             gp=gpar(lwd = 1, col = "darkred")) # lwd = line width, col = line color
  
  ## This code adds the x and y axis lines
  # at = is a numeric vector with the x-axis locations for tick marks
  grid.xaxis(at = 
               c(seq(1, flankBP, length.out = 3),
                 flankBP + seq(1, motifWidth),
                 flankBP + motifWidth + seq(1, flankBP, length.out = 3)),
             label = c(-(flankBP + 1 - seq(1, flankBP + 1, length.out = 3)),
                       rep("", motifWidth),
                       seq(0, flankBP, len = 3)))
  grid.yaxis()
  
  ## Adds the dashed line across the x-axis horizontally (motif hashes)
  grid.lines(x=c(flankBP, flankBP, 0), y=c(0, max(insVector), ylim[2]),
             default.units="native", gp=gpar(lty=2))
  
  ##
  grid.lines(x=c(flankBP + motifWidth + 1, flankBP + motifWidth + 1, totalBP),
             y=c(0, max(insVector), ylim[2]),
             default.units="native", gp=gpar(lty=2))
  upViewport()
  
  ##
  vp2 <- viewport(y=.9, height=.2,
                  xscale=c(0, totalBP + 1),
                  name="motif")
  pushViewport(vp2)
  upViewport()
  
  ##
  legvp <- viewport(x=0.5,
                    y=0.5,
                    width=convertX(unit(1, "lines"), unitTo="npc"),
                    height=convertY(unit(1, "lines"), unitTo="npc"),
                    just=c("right", "top"), name="legendWraper")
  pushViewport(legvp)
  upViewport()
  
  ##
  grid.text(plotTitle,
            y=unit(1, "npc")-convertY(unit(1, "lines"), unitTo="npc"),
            gp=gpar(cex=1.2, fontface="bold"))
  upViewport()
  
  ## Add the x and y axis labels to the image
  grid.text(xlab, y=unit(1, 'lines'))
  grid.text(ylab, x=unit(1, 'line'), rot = 90)
  
} # end plotInsProb function

## Make graph of the raw peak sites
#svgPath <- paste0(dirPath, "footprints/graphs/peaks/", sampleName, ".", geneName, ".", "motif", a, ".rawpeak.sites.svg")
#svg(file = svgPath)
#cat("Saving peaks footprint image at path:", svgPath, "\n")
#plotTitle <- paste0(sampleName, ".", geneName, ".", "motif", a, ".rawpeaks")
#plotInsProb(plotTitle = plotTitle, motifWidth = motifWidth, motifPWM = PWM, insVector = insVector)
#dev.off()