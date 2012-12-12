rm(list=ls(all=T))
library("shiny")
library("brew")
library("flowCore")
library("Rclusterpp")


dataDir = "/Users/rbruggner/Desktop/notime/citrusTestRun/"
sapply(list.files("/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI/guiFunctions/",pattern=".R",full.names=T),source)
source("/Users/rbruggner/Desktop/work/citrus/R/citrus.launchUI.R")

fileList = list.files(dataDir,pattern=".fcs",ignore.case=T)
fileGroupAssignments = rep("",length(fileList))  
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)


result = tryCatch({
  runApp(appDir="/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI",launch.browser=T)
}, warning = function(w) {
   cat(as.character(w))
}, error = function(e) {
    stop(paste("Unexpected Error:",e))
}, finally = {
  
})

