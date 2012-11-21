citrus.readFCSSet = function(dataDirectory,conditionFileList,conditionSampleSize=NULL,transformColumns=NULL){
  data = list();
  fileCounter = 1;
  conditions = colnames(conditionFileList)
  for (i in 1:length(conditions)){
    cat(paste("Reading Condition ",conditions[i],"\n"));
    conditionData = list();
    for (fileName in conditionFileList[,i]){
      filePath = paste(dataDirectory,fileName,sep="");
      if (!file.exists(filePath)){
        stop(paste("File",filePath,"not found."));
      }
      cat(paste("\tReading file ",fileName,"\n"));
      fcsData = exprs(read.FCS(filePath));
      fcsData = cbind(fcsData,fileEventNumber=1:nrow(fcsData),fileId=fileCounter);
      fileCounter=fileCounter+1;
      if (!is.null(transformColumns)){
        fcsData[,transformColumns] = asinh(fcsData[,transformColumns]/5);
      }
      if ((!is.null(conditionSampleSize))&&(conditionSampleSize[i]<nrow(fcsData))){
        cat(paste("\tSampling",conditionSampleSize[i],"events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),conditionSampleSize[i])),] 
      }
      conditionData[[fileName]] = fcsData
    }
    data[[conditions[i]]] = do.call("rbind",conditionData)
    rm(conditionData)
    gc();
  }
  
  results = list(data=do.call("rbind",data),fileIds=matrix(1:(fileCounter-1),ncol=ncol(conditionFileList),dimnames=list(c(),conditions)),fileNames=fileNames)
  class(results) = "citrusDataObject"
  # HOW DO I MAKE IT PRINT OUT SUMMARIES? 
  return(results);
}

citrus.indexof = function(x,y){
  which(y==x)
}

citrus.leavoutFold = function(x,y,leaveoutSize){
  groupCounts = table(y)/length(y)*leaveoutSize
  leaveout = list()
  for (group in unique(y)){
    leaveout[[group]] = sample(which(y==group),groupCounts[which(names(groupCounts)==group)])
  }
  return(as.vector(unlist(leaveout)))
}