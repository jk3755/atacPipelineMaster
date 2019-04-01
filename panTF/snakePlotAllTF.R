
## Load Packages
loadLibrary <- function(lib){
  for( i in lib ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )}}}
##
loadLibrary(c("ggplot2", "ggsci"))

## Install libraries, if necessary
install.packages("ggplot2")
install.packages("ggsci")
##
library(ggplot2)
library(ggsci)



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

##
ggplot(data, aes(footprint, flank, color=background)) + 
  geom_point() + 
  scale_fill_gsea()



## Plot 2 - Footprint depth against flanking accessibility, for all TFs, dots colored by VIPER protein activity level
