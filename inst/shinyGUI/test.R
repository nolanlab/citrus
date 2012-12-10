library("shiny")
library("flowCore")
options(shiny.trace=FALSE)

guiPath = "/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI/"
source(paste(guiPath,"initFunctions.R",sep=""))

dataDir = "/Users/rbruggner/Desktop/work/citrus/data/syntheticData/train/unstim/"
fileList = list.files(dataDir)
fileGroupAssignments = rep("",length(fileList))  
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)


runApp("/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI")