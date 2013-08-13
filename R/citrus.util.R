# FILE SAMPLE SIZE SHOULD BE A NAMED VECTOR OR LIST OR SOMETHING THAT'S EASY TO EXTRACT BY NAME
citrus.readFCSSet = function(dataDir,fileList,conditions,fileSampleSize=NULL,transformCols=NULL,transformFactor=5,emptyValue=T){
  data = list();
  fileCounter = 1;
  fileNames = c();
  fileChannelNames = list();
  fileReagentNames = list();
  for (i in 1:length(conditions)){
    cat(paste("Reading Condition ",conditions[i],"\n"));
    conditionData = list();
    fileChannelNames[[conditions[i]]] = list();
    fileReagentNames[[conditions[i]]] = list();
    for (fileName in fileList[,conditions[i]]){
      fileNames[fileCounter]=fileName
      filePath = file.path(dataDir,fileName);
      if (!file.exists(filePath)){
        stop(paste("File",filePath,"not found."));
      }
      cat(paste("\tReading file ",fileName,"\n"));
      suppressWarnings((fcsFile = read.FCS(filePath,emptyValue=emptyValue)))
      fcsData=exprs(fcsFile)
      fileChannelNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$name)
      fileReagentNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$desc)
      fcsData = cbind(fcsData,fileEventNumber=1:nrow(fcsData),fileId=fileCounter);
      fileCounter=fileCounter+1;
      if (!is.null(transformCols)){
        fcsData[,transformCols] = asinh(fcsData[,transformCols]/transformFactor);
      }
      if ((!is.null(fileSampleSize))&&(fileSampleSize<nrow(fcsData))){
        cat(paste("\tSampling",fileSampleSize,"events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),fileSampleSize)),] 
      }
      conditionData[[fileName]] = fcsData
    }
    data[[conditions[i]]] = do.call("rbind",conditionData)
    rm(conditionData)
    gc();
  }
  
  results = list(data=do.call("rbind",data),fileIds=matrix(1:(fileCounter-1),ncol=length(conditions),dimnames=list(c(),conditions)),fileNames=fileNames,fileChannelNames=fileChannelNames,fileReagentNames=fileReagentNames)
  class(results) = "citrusDataObject"
  # HOW DO I MAKE IT PRINT OUT SUMMARIES? 
  return(results);
}

citrus.assembleHandGates = function(dataDir,filePopulationList,conditionComparaMatrix=NULL,fileSampleSize=NULL,transformColumns=NULL,transformFactor=5){
  
  if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
    # Perform Internal Consistency Checks!
  } else {
    allConditions = names(filePopulationList)
  }
  
  results = list();

  for (conditions in allConditions){
    
    cda = data.frame();
    eventCounter=0;
    populationCounter=1;  
    fileNames=list()  
    clusterAssignments = list();
    fileId=1;  
    fileChannelNames = list()
    fileReagentNames = list()
    fileIds = list();
    populationClusterMap = list();
    for (condition in conditions){
      
      nPatients = nrow(filePopulationList[[condition]])
      baseFileId = (match(condition,conditions)-1)*nPatients
      fileIds[[condition]] = (baseFileId+1):(baseFileId+nPatients)
      for (populationName in colnames(filePopulationList[[condition]])){
        populationEvents = c();
        for (patientId in 1:nPatients){
          fcsFileName = filePopulationList[[condition]][patientId,populationName]
          fileId = baseFileId+patientId
          print(paste("Reading",fcsFileName))
          fcsFile = read.FCS(file.path(dataDir,fcsFileName))
          fileChannelNames[[fcsFileName]] = colnames(fcsFile)
          fileReagentNames[[fcsFileName]] = pData(parameters(fcsFile))$desc
          fcsFileData = exprs(fcsFile)
          if (!is.null(fileSampleSize)){
            if (fileSampleSize < nrow(fcsFileData)){
              fcsFileData = fcsFileData[sample(1:nrow(fcsFileData),fileSampleSize),]
            }
          }
          lastEvent = nrow(cda)
          absoluteEventIds = (lastEvent+1):(eventCounter+nrow(fcsFileData))
          cda = rbind(cda,data.frame(fcsFileData,fileId=fileId))
          populationEvents = c(populationEvents,absoluteEventIds)
          eventCounter = eventCounter+nrow(fcsFileData)
          fileNames[[fileId]]=fcsFileName
          fileId=fileId+1;
        }
        
        
        if (populationName %in% names(populationClusterMap)){
          clusterAssignments[[populationClusterMap[[populationName]]]] = c(clusterAssignments[[populationClusterMap[[populationName]]]],populationEvents)
        } else {
          clusterAssignments[[populationCounter]]=populationEvents;
          populationClusterMap[[populationName]]=populationCounter
          populationCounter = populationCounter+1;
        }
        
      }
    }
    
    if (!is.null(transformColumns)){
      cda[,transformColumns] = cda[,transformColumns]/transformFactor
    }
    
    results[[paste(conditions,collapse="_vs_")]] = list(folds="all",
                                                        foldsClusterAssignments=list(all=clusterAssignments),
                                                        conditions=conditions,
                                                        citrus.dataArray=list(data=as.matrix(cda),
                                                                              fileNames=unlist(fileNames,use.names=F),
                                                                              fileIds=do.call("cbind",fileIds),
                                                                              fileChannelNames=fileChannelNames,
                                                                              fileReagentNames=fileReagentNames),
                                                        clusterPopulationNames=names(populationClusterMap));
  }
  return(results)
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

citrus.formatDecimal = function(x){
  sprintf("%1.2f", x)
}

citrus.convertConditionMatrix = function(conditionMatrix){
  conditions = list();
  for (i in 1:nrow(conditionMatrix)){
    for (j in 1:ncol(conditionMatrix)){
      if (conditionMatrix[i,j]){
        conditions = append(conditions,list(unique(c(rownames(conditionMatrix)[i],colnames(conditionMatrix)[j]))))
      }
    }
  }
  return(conditions)
}

citrus.version = function(){
  return("0.05")
}

citrus.fileEventCount = function(dataDir,emptyValue=T){
  lengths = list();
  for (fcsFile in list.files(dataDir,pattern=".fcs",ignore.case=T)){
    print(paste("Reading",fcsFile))
    lengths[[fcsFile]] = suppressWarnings(dim(read.FCS(file.path(dataDir,fcsFile),emptyValue=emptyValue)))
  }
  return(do.call("rbind",lengths))
}

citrus.featureTypes = function(){
  return(c("densities","medians","emDists"))
}
citrus.familyList = function(){
  return(c("twoClass","survival","quantiative"))
}
citrus.modelTypes = function(){
  return(c("pamr","glmnet","sam"))
}
