runCitrus = FALSE;

preload=F
# Choose any file from the appropriate directory

if (!exists("dataDir")){
  dataDirFile = file.choose()
  if (is.null(dataDirFile)){
    stop("File Selection Canceled")
  }
  dataDir = dirname(dataDirFile)
} else {
  if (!file.exists(dataDir)){
    stop(paste("Directory",dataDir,"not found. Exiting."))
  }
  dataDirFile = dataDir;
}

if (basename(dataDirFile)=="citruskey.csv"){
  tryCatch({
    keyFile = read.csv(dataDirFile,header=T)
  }, warning = function(w) {
    stop(simpleWarning(paste("Error Reading input file:",w)))
  }, error = function(e) {
    stop(simpleError(paste("Error Reading input file:",e)))
  }, finally = {
    if (!("labels" %in% colnames(keyFile))){
      stop("Error reading citrus key: 'labels' column not found");
    } else {
      labelCol = which(colnames(keyFile)=="labels")
    }
    preload=T
    dataDir = dirname(dataDirFile)
  });
}
    
cat(paste("Launching citrus interface with target directory:",dataDir,"\n"));

# Lousy Fix to keep open mp processes from dying in mac os x
#if (Sys.info()[["sysname"]]=="Darwin"){
#  cat("Mac OS X Detected. Setting number of OpenMP threads to 1.\n")
#  Rclusterpp.setThreads(threads=1)
#  Sys.setenv("OMP_NUM_THREADS"=1)
#}

# Comment to True to debug
options(shiny.trace=F)

# Get directory 

if (!preload){
  # Get list of sample files
  fileList = list.files(file.path(dataDir),pattern=".fcs",ignore.case=T)
  # This should get fixed...
  fileGroupAssignments = rep("",length(fileList))
  # Pre-read list of columns measured in each file
} else {
  fileList = as.vector(unlist(keyFile[,-labelCol]))
  fileGroupAssignments = as.vector(rep(keyFile[,labelCol],ncol(keyFile[,-labelCol])))
} 

fileCols = lapply(fileList,citrus.getFileCols,dataDir=dataDir)


disableInput <- function(x) {
  if (inherits(x, 'shiny.tag')) {
    if (x$name %in% c('input', 'select'))
      x$attribs$disabled <- 'disabled'
    x$children <- disableInput(x$children)
  }
  else if (is.list(x) && length(x) > 0) {
    for (i in 1:length(x))
      x[[i]] <- disableInput(x[[i]])
  }
  x
}