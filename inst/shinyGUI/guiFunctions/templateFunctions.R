
writeRunCitrusFile = function(input,templateFile=NULL){
  templateData = as.list(input)
  templateData[["dataDir"]]=dataDir
  outputDir = paste(dataDir,"citrusOutput",sep="")
  if (!file.exists(outputDir)){
    dir.create(file.path(dataDir,"citrusOutput"),showWarnings=F)
  }
  brew(
    file="/Users/rbruggner/Desktop/work/citrus/inst/shinyGUI/guiFunctions/runCitrus.template",
    output=paste(outputDir,"runCitrus.R",sep=.Platform$file.sep)
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