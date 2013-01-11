citrus.launchUI = function(dataDirectory=NULL){  
  
  library("shiny")
  library("brew")
  
  if (!is.null(dataDirectory)){
    dataDir <<-dataDirectory
  }
  
  #sapply(list.files(file.path(system.file(package = "citrus"),"shinyGUI","guiFunctions"),pattern=".R",full.names=T),source)
  
  res = tryCatch({
    runApp(appDir=file.path(system.file(package = "citrus"),"shinyGUI"),launch.browser=T)  
  }, warning = function(w){
    cat(paste(w,"\n"));
  },error = function(e){
    stop(print(paste("Unexpected Error:",e)))
  }, finally = {
    
  })
  
  if (runCitrus){
    outputPath = file.path(dataDir,"citrusOutput")
    setwd(outputPath)
    runFile = file.path(outputPath,"runCitrus.R")
    cat(paste("Running Citrus File:",runFile,"\n"))
    logFilePath = file.path(outputPath,"citrusOutput.log")
    cat(paste("Logging output to:",logFilePath,"\n"))
    logFile = .logOn(logFilePath)
    
    source(runFile)
    
    .logOff(logFile=logFile)
  }
  return(paste("Citrus Output in:",outputPath))  
}

citrus.getClusterCols = function(fileName,dataDir){
  fcsFile = suppressWarnings(read.FCS(file.path(dataDir,fileName),which.lines=1))
  return(flowCore::colnames(fcsFile))
}

.logOn = function(filePath,messages=F){
  logFile = file(filePath,open = "wt")
  sink(logFile,split=T)
  if (messages){
    sink(logFile,type="message")
  }
  return(logFile)
}

.logOff = function(logFile,messages=F){
  if (messages){
    sink(type="messsage")
  }
  sink()
  close(logFile)
}