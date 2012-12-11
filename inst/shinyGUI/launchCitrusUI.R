###########################################################
# This script launches the shiny web UI for Citrus
# To launch UI from a command line,type:
# >R -f /path/to/launchCitrusUI.R /path/to/fcs/files/
###########################################################
rm(list=ls())
LIBRARY_PATH=NULL
library("citrus",lib.loc=LIBRARY_PATH)
library("shiny",lib.locLIBRARY_PATH)
library("brew")
library("snow")

# Load functions to retreive initial variables
initFunctionPath = file.path(system.file(package = "citrus"),"shinyGUI","guiFunctions","launcherInitFunctions.R")
source(initFunctionPath)

# Comment to True to debug
options(shiny.trace=FALSE)

# Get directory 
dataDir = "/Users/rbruggner/Desktop/notime/citrusTestRun"
#dataDir = commandArgs()[4]
if (!file.exists(dataDir)){
  stop(paste("Directory",dataDir,"not found. Exiting."))
}

# Get list of sample files
fileList = list.files(file.path(dataDir),pattern=".fcs",ignore.case=T)

# This should get fixed...
fileGroupAssignments = rep("",length(fileList))  

# Pre-read list of columns measured in each file
fileCols = lapply(fileList,getClusterCols,dataDir=dataDir)

# Launch Shiny UI in SNOW thread so we can quit it with out shutting down all of R.
guiThread = makeCluster(1,type="SOCK")
clusterEvalQ(guiThread,library("brew"))
clusterExport(guiThread,c("dataDir","fileList","fileGroupAssignments","fileCols"))
clusterCall(guiThread,runApp,appDir=file.path(system.file(package = "citrus"),"shinyGUI"),launch.browser=T)

# Call newly written citrus file
source(file.path(dataDir,"citrusOutput","runCitrus.R"))