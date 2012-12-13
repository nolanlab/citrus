citrus.launchUI = function(dataDirectory=NULL){  
  
  library("shiny")
  library("brew")
  
  if (!is.null(dataDirectory)){
    dataDir <<-dataDirectory
  }
  
  sapply(list.files(file.path(system.file(package = "citrus"),"shinyGUI","guiFunctions"),pattern=".R",full.names=T),source)
  
  res = tryCatch({
    runApp(appDir=file.path(system.file(package = "citrus"),"shinyGUI"),launch.browser=T)  
  }, warning = function(w){
    cat(paste(w,"\n"));
  },error = function(e){
    stop(print(paste("Unexpected Error:",e)))
  }, finally = {
    
  })
  runFile = file.path(dataDir,"citrusOutput","runCitrus.R")
  cat(paste("Running Citrus File:",runFile,"\n"))
  source(runFile)
  return(paste("Citrus Output in:",file.path(dataDir,"citrusOutput")))
}

citrus.getClusterCols = function(fileName,dataDir){
  fcsFile = suppressWarnings(read.FCS(file.path(dataDir,fileName),which.lines=1))
  return(flowCore::colnames(fcsFile))
}