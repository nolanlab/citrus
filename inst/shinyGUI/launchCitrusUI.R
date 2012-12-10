###########################################################
# This script launches the shiny web UI for Citrus
# To launch UI from a command line,type:
# >R -f /path/to/launchCitrusUI.R /path/to/fcs/files/
###########################################################

rm(list=ls())
LIBRARY_PATH=NULL
library("citrus",lib.loc=LIBRARY_PATH)
library("shiny",lib.locLIBRARY_PATH)

# Load functions to retreive initial variables
initFunctionPath = system.file(paste("shinyGUI","initFunctions.R",sep=.Platform$file.sep), package = "citrus")
source(initFunctionPath)

# Comment to True to debug
options(shiny.trace=FALSE)

# Get directory 
dataDir = commandArgs()[4]
if (!file.exists(dataDir)){
  stop(paste("Directory",dataDir,"not found. Exiting."))
}

# Get list of sample files
fileList = list.files(dataDir)

# This should get fixed...
fileGroupAssignments = rep("",length(fileList))  

# Pre-read list of columns measured in each file
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)

# Launch Shiny UI
runApp(system.file(paste("shinyGUI",sep=.Platform$file.sep), package = "citrus"),launch.browser=T)
