
# Choose any file from the appropriate directory
if (!exists("dataDir")){
  dataDirFile = file.choose()
  if (is.null(dataDirFile)){
    stop("File Selection Canceled")
  } else {
    dataDir = dirname(dataDirFile)
  }  
}

cat(paste("Launching citrus interface with target directory:",dataDir,"\n"));

# Lousy Fix to keep open mp processes from dying in mac os x
if (Sys.info()[["sysname"]]=="Darwin"){
  cat("Mac OS X Detected. Setting number of OpenMP threads to 1.\n")
  Rclusterpp.setThreads(threads=1)  
}

# Comment to True to debug
options(shiny.trace=FALSE)

# Get directory 
if (!file.exists(dataDir)){
  stop(paste("Directory",dataDir,"not found. Exiting."))
}

# Get list of sample files
fileList = list.files(file.path(dataDir),pattern=".fcs",ignore.case=T)

# This should get fixed...
fileGroupAssignments = rep("",length(fileList))

# Pre-read list of columns measured in each file
fileCols = lapply(fileList,citrus.getClusterCols,dataDir=dataDir)

