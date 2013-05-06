citrus.foldCluster = function(foldMembership,citrus.dataArray,clusterCols,conditions){
  if ((length(foldMembership)==1) && (foldMembership=="all")){
    cat("Clustering All\n")
    includeFileIds = as.vector(citrus.dataArray$fileIds[,conditions])
  } else {
    includeFileIds = as.vector(citrus.dataArray$fileIds[-foldMembership,conditions])
    cat(paste("Clustering fileIds",paste(includeFileIds,collapse=", "),"\n"))
  }
  return(citrus.cluster(citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% includeFileIds,clusterCols]))
}

citrus.cluster = function(data){
  cat(paste("Clustering",nrow(data),"events\n"));
  return(Rclusterpp.hclust(data))
}

citrus.calculateCompleteHierarchicalMembership = function(clustering){
    mergeOrder = clustering$merge
    return(lapply(as.list(1:nrow(mergeOrder)),citrus.traverseMergeOrder,mergeOrder=mergeOrder))
}

citrus.mapFoldDataToClusterSpace = function(index,citrus.dataArray,foldClusterAssignments,folds,conditions,clusterCols){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    return(NULL)
  }
  foldFileIds = as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
  leftoutFileIds = as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
  return(citrus.mapDataToClusterSpace(data=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%foldFileIds,clusterCols],clusterAssignments=foldClusterAssignments[[index]],newData=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%leftoutFileIds,clusterCols]))  
}

citrus.mapDataToClusterSpace = function(data,clusterAssignments,newData){  
  nnMap = citrus.assignToCluster(tbl=newData,cluster_data=data,cluster_assign=rep(1,nrow(data)))
  return(lapply(clusterAssignments,citrus.mapNeighborsToCluster,nearestNeighborMap=nnMap))
}  

citrus.mapNeighborsToCluster = function(clusterMembers,nearestNeighborMap){
  which(nearestNeighborMap %in% clusterMembers)
}

citrus.traverseMergeOrder = function(node,mergeOrder){
  left = mergeOrder[node,1]
  right = mergeOrder[node,2]
  if (left<0){
    leftAtom = -left;  
  } else {
    leftAtom = citrus.traverseMergeOrder(left,mergeOrder);
  }
  
  if (right<0){
    rightAtom = -right;
  } else {     
    rightAtom = citrus.traverseMergeOrder(right,mergeOrder); 
  }
  return(c(leftAtom,rightAtom));
}


citrus.preCluster = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,transformCols=NULL,conditionComparaMatrix=NULL,balanceFactor=NULL,transformFactor=5,...){
  
  addtlArgs = list(...)
  
  res=list()
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  emptyValue=T
  if ("emptyValue" %in% names(addtlArgs)){
    emptyValue = addtlArgs[["emptyValue"]]
  }

  if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
  } else {
    allConditions = as.list(colnames(fileList))
  }
  
  if (is.null(balanceFactor)){
    balanceFactor = as.factor(sample(c(0,1),nrow(fileList),replace=T))
  } else if (levels(balanceFactor)==1){
    balanceFactor = as.factor(sample(c(0,1),length(balanceFactor),replace=T))
  }
  if (nFolds=="all"){
    folds = list()
    nAllFolds=1
  } else {
    folds = pamr:::balanced.folds(y=balanceFactor,nfolds=nFolds)
    nAllFolds = nFolds+1  
  }
  folds[[nAllFolds]]="all"
  
  for (conditions in allConditions){
    cat(paste("Clustering Condition",paste(conditions,collapse=" vs "),"\n"))
    
    citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize,transformFactor=transformFactor,emptyValue=emptyValue)
    
    foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions)
    cat("Assigning Events to Clusters\n")
    if ("snowCluster" %in% names(addtlArgs)){
      foldsClusterAssignments = parLapply(addtlArgs[["snowCluster"]],foldsCluster,citrus.calculateCompleteHierarchicalMembership)
    } else {
      foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)  
    }
    
    resObj = list(folds=folds,foldsCluster=foldsCluster,foldsClusterAssignments=foldsClusterAssignments,conditions=conditions,citrus.dataArray=citrus.dataArray,clusterColumns=clusterCols)
    if (nFolds!="all"){
      cat("Assigning Leftout Events to Clusters\n")
      if ("snowCluster" %in% names(addtlArgs)){
        leftoutClusterAssignments = parLapply(addtlArgs[["snowCluster"]],1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols)
      } else {
        leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols)  
      }
      
      resObj$leftoutClusterAssignments=leftoutClusterAssignments
    }
    save(resObj,file=file.path(outputDir,paste("citrus.Cluster.",paste(conditions,collapse="_vs_"),".rDat",sep="")))
    res[[paste(conditions,collapse="_vs_")]]=resObj
  }
  return(res)
}
