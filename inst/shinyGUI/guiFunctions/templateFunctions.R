
writeRunCitrusFile = function(input,templateFile=NULL){
  templateData = as.list(input)
  templateData[["dataDir"]]=dataDir
  outputDir = file.path(dataDir,"citrusOutput")
  if (!file.exists(outputDir)){
    dir.create(file.path(dataDir,"citrusOutput"),showWarnings=F)
  }
  brew(
    file=file.path(system.file(package="citrus"),"shinyGUI","guiFunctions","runCitrus.template"),
    output=file.path(outputDir,"runCitrus.R")
    )
}

convertToLabeledFileList = function(x){
    df = data.frame(labels=c(),defaultCondition=c())
    for (group in names(x)){
      df = rbind(df,data.frame(labels=rep(group,length(x[[group]])),defaultCondition=x[[group]]))
    }
    df[,"labels"]=as.factor(df[,"labels"])
    return(df)
}

convertColToDefinition = function(colname,df){
  paste(colname,"=c(",paste(sapply(df[,colname],stringQuote),collapse=","),")",sep="")
}