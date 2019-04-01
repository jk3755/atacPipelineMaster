
library(ggplot2)

## Plot 1 - Footprint depth against flanking accessibility, for all TFs, dots colored by RNA expression level


## List all the files that will be analyzed
fileList <- list.files("C:\\Users\\jsk33\\Desktop\\fpdata\\", full.names = TRUE)
numFiles <- length(fileList)

##
dataMatrix <- matrix(data = NA, nrow = numFiles, ncol = 3)
colnames(dataMatrix) <- c("Flank", "Motif", "Background")

##
for (a in 1:numFiles){
  
  load(fileList[a])
  temp <- footprintData[["motif1"]][["rawFootprintMetrics"]]
  ##
  dataMatrix[a,1] <- mean(temp[,4])
  dataMatrix[a,2] <- mean(temp[,5])
  dataMatrix[a,3] <- mean(temp[,1])
  
}

##
data <- data.frame(flank = dataMatrix[,1], footprint = dataMatrix[,2], background = dataMatrix[,3])

##
ggplot(data, aes(footprint, flank, color=background)) + 
  geom_point() + 
  scale_color_gradient(low="blue", high="red")





## Plot 2 - Footprint depth against flanking accessibility, for all TFs, dots colored by VIPER protein activity level
