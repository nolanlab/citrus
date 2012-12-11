rm(list=ls(all=T))
library("flowCore")
library("shiny")
library("brew")

options(shiny.trace=FALSE)

guiPath = "/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI/"
sapply(list.files(paste(guiPath,"guiFunctions",sep=""),pattern=".R",full.names=T),source)

dataDir = "/Users/rbruggner/Desktop/notime/citrusTestRun/"
fileList = list.files(dataDir)
fileGroupAssignments = rep("",length(fileList))  
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)

runApp("/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI",launch.browser=T)

#library("snow")
#thread = makeCluster(1,type="SOCK")
#clusterEvalQ(thread,library("shiny"))
#clusterEvalQ(thread,)

