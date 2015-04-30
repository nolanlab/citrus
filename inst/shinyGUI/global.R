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
    keyFile = read.csv(dataDirFile,header=T,stringsAsFactors=F)
  }, warning = function(w) {
    stop(simpleWarning(paste("Error Reading input file:",w)))
  }, error = function(e) {
    stop(simpleError(paste("Error Reading input file:",e)))
  }, finally = {
    if ("class" %in% colnames(keyFile)){
      labelCol = which(colnames(keyFile)=="class")
      family="classification"
    } else if ("endpoint_value" %in% colnames(keyFile)){
      labelCol = which(colnames(keyFile)=="endpoint_value")
      family="continuous"
    } else {
      stop("Error reading citrus key: endpoint column 'class' or 'endpoint_value' not found");
    }
    preload=T
    dataDir = dirname(dataDirFile)
  });
}
    
cat(paste("Launching citrus interface with target directory:",dataDir,"\n"));


# Comment to True to debug
options(shiny.trace=F)

# Get directory 
if (!preload){
  # Get list of sample files
  fileList = list.files(file.path(dataDir),pattern=".fcs",ignore.case=T)
  conditionFiles = data.frame(defaultCondition=fileList)
  
  if (length(fileList)==0){
    stop(paste0("\nNo FCS files found in  ",dataDir,". Please ensure files have a '.fcs' or '.FCS' extension."))
  }
  
  # This should get fixed...
  fileGroupAssignments = rep("",length(fileList))
  # Pre-read list of columns measured in each file
} else {
  fileList = as.vector(unlist(keyFile[,-labelCol]))
  conditionFiles = keyFile[,-labelCol,drop=F]
  if (family=="classification"){
    fileGroupAssignments = as.vector(rep(keyFile[,labelCol],ncol(keyFile[,-labelCol])))  
  } else if (family=="continuous") {
    sampleEndpointValues = keyFile[,labelCol]
  } else {
    stop(paste("Unknown family:",family))
  }
  
} 

cat("\nScanning parameters in FCS files\n")
fileCols = lapply(fileList,citrus.getFileParameters,dataDir=dataDir)
fileColLength = sapply(fileCols,length)
cat("\nNumber of parameters per file:\n")
cat(paste0(fileList,": ",fileColLength,"\n"))
if (length(unique(fileColLength))>1){
  stop("\nAll FCS files must have the same number of channels.\n")
}
  
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
