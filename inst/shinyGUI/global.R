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
}

if (basename(dataDirFile)=="citruskey.csv"){
  # Preload file data and labels
  preload=T
  keyFile = tryCatch({
    read.csv(dataDirFile,header=T,stringsAsFactors=F)
  }, warning = function(w) {
    stop(simpleWarning(paste("Error Reading input file:",w)))
  }, error = function(e) {
    stop(simpleError(paste("Error Reading input file:",e)))
  });
  
  if ("class" %in% colnames(keyFile)){
    labelCol = which(colnames(keyFile)=="class")
    family="classification"
    fileGroupAssignments = as.vector(rep(keyFile[,labelCol],ncol(keyFile[,-labelCol])))  
  } else if ("endpoint_value" %in% colnames(keyFile)){
    labelCol = which(colnames(keyFile)=="endpoint_value")
    family="continuous"
  } else {
    stop("Error reading citrus key: endpoint column 'class' or 'endpoint_value' not found");
  }
  
  
  # Set directory where data live
  dataDir = dirname(dataDirFile)
  
  # Set labels 
  labels = keyFile[,labelCol]
  
  # Get list of initial conditions that are available in file
  conditions = colnames(keyFile[,-labelCol])
  
  # List of all files, regardless of condition
  fileList = as.vector(unlist(keyFile[,-labelCol]))
  
  # Files by condition
  conditionFiles = keyFile[,-labelCol,drop=F]
  
} else {
  # No preload. Assume family classification until new dynamic interface
  family="classification"
  
  # List of all files, regardless of condition
  fileList = list.files(file.path(dataDir),pattern=".fcs",ignore.case=T)
  if (length(fileList)==0){
    stop(paste0("\nNo FCS files found in  ",dataDir,". Please ensure files have a '.fcs' or '.FCS' extension."))
  }
  
  # Assign all files to default condition
  conditionFiles = data.frame(defaultCondition=fileList)
  conditions = "defaultCondition"
  
  # Ugly hack.
  # Assign all files to the "" group until assigned to a specific group name
  # This should get fixed...
  fileGroupAssignments = rep("",length(fileList))
}


# Begin launch sequence
# Emit basic info and run consistency checks
options(shiny.trace=F)
cat(paste("Regression family:",family,"\n"))
cat(paste("Launching citrus interface with target directory:",dataDir,"\n"));

# Perform basic parameter checks
cat("\nScanning parameters in FCS files\n")
fileCols = lapply(fileList,citrus.getFileParameters,dataDir=dataDir)
fileCheckResult = citrus.checkFileParameterConsistency(dataDir,fileList)
if (fileCheckResult!=0){
  stop("\nInconsistent parameters in FCS files.\n")
}
  
