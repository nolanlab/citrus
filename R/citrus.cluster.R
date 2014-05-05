citrus.mapFileDataToClustering = function(dataDir,newFileList,preClusterResult,fileSampleSize,mappingColumns=NULL,transformCols=NULL,transformCofactor=5,...){
  conditions = strsplit(names(preClusterResult),"_vs_")
  mappingResult = list()
  for (conditionName in names(preClusterResult)){
    conditions = strsplit(conditionName,"_vs_")[[1]]
    mappingResult[[conditionName]]$citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=newFileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize,transformCofactor=transformCofactor,...)
    
    if (is.null(mappingColumns)){
      mappingColumns = preClusterResult[[conditionName]]$clusterColumns
    }
    
    cat(paste("Mapping",nrow(mappingResult[[conditionName]]$citrus.dataArray$data),"events to existing cluster space.\n"));
    mappingResult[[conditionName]]$foldsClusterAssignments = list()
    mappingResult[[conditionName]]$foldsClusterAssignments[["mappingResult"]] = citrus.mapDataToClusterSpace(data=preClusterResult[[conditionName]]$citrus.dataArray$data[,mappingColumns],
                                 clusterAssignments=preClusterResult[[conditionName]]$foldsClusterAssignments[[which(preClusterResult[[conditionName]]$folds=="all")]],
                                 newData=mappingResult[[conditionName]]$citrus.dataArray$data[,mappingColumns],...)
    mappingResult[[conditionName]]$folds = list("mappingResult")
    mappingResult[[conditionName]]$conditions = conditions
    mappingResult[[conditionName]]$clusterColumns = preClusterResult[[conditionName]]$clusterColumns
  }
  return(mappingResult)
}


citrus.foldCluster = function(foldMembership,citrus.dataArray,clusterCols,conditions,...){
  if ((length(foldMembership)==1) && (foldMembership=="all")){
    cat("Clustering All\n")
    includeFileIds = as.vector(citrus.dataArray$fileIds[,conditions])
  } else {
    includeFileIds = as.vector(citrus.dataArray$fileIds[-foldMembership,conditions])
    cat(paste("Clustering fileIds",paste(includeFileIds,collapse=", "),"\n"))
  }
  return(citrus.cluster(citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% includeFileIds,clusterCols],...))
}

citrus.cluster = function(data,...){
  cat(paste("Clustering",nrow(data),"events\n"));
  
  clusterType="hierarchical"
  if ("clusterType" %in% names(list(...))){
    clusterType = list(...)[["clusterType"]]
  }
  if (clusterType=="hierarchical"){
    return(Rclusterpp.hclust(data))  
  } else {
    stop(paste("Don't know how to cluster using",clusterType))
  }
  
}

citrus.calculateCompleteHierarchicalMembership = function(clustering,...){
  clusterType="hierarchical"
  if ("clusterType" %in% names(list(...))){
    clusterType = list(...)[["clusterType"]]
  }
  if (clusterType=="hierarchical"){
    mergeOrder = clustering$merge
    if ("mc.cores" %in% names(list(...))){
      return(mclapply(as.list(1:nrow(mergeOrder)),citrus.traverseMergeOrder,mergeOrder=mergeOrder,mc.cores=list(...)[["mc.cores"]]))
    } else {
      return(lapply(as.list(1:nrow(mergeOrder)),citrus.traverseMergeOrder,mergeOrder=mergeOrder))
    }
  } else {
    stop(paste("Don't know how to caclulate complete clustering for clusterType",clusterType))
  }
}

citrus.mapFoldDataToClusterSpace = function(index,citrus.dataArray,foldClusterAssignments,folds,conditions,clusterCols,...){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    return(NULL)
  }
  foldFileIds = as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
  leftoutFileIds = as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
  return(citrus.mapDataToClusterSpace(data=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%foldFileIds,clusterCols],clusterAssignments=foldClusterAssignments[[index]],newData=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%leftoutFileIds,clusterCols],...))  
}

citrus.mapDataToClusterSpace = function(data,clusterAssignments,newData,...){  
  nnMap = citrus.assignToCluster(tbl=newData,cluster_data=data,cluster_assign=rep(1,nrow(data)))
  addtlArgs = list(...)
  if ("mc.cores" %in% names(addtlArgs)){
    return(mclapply(clusterAssignments,citrus.mapNeighborsToCluster,nearestNeighborMap=nnMap,mc.cores=addtlArgs[["mc.cores"]]))
  } else {
    return(lapply(clusterAssignments,citrus.mapNeighborsToCluster,nearestNeighborMap=nnMap))
  }
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

citrus.getClusterDecendants = function(node,mergeOrder){
  left = mergeOrder[node,1]
  right = mergeOrder[node,2]
  if (left>0){
    left = c(left,citrus.getClusterDecendants(left,mergeOrder))
  } else {
    left = c()
  }
  if (right>0){
    right = c(right,citrus.getClusterDecendants(right,mergeOrder))
  } else {
    right = c()
  }
  return(c(left,right))
}

citrus.getClusterAncestors = function(node,mergeOrder){
  parent = which(mergeOrder==node,arr.ind=T)[1]
  if (is.na(parent)){
    return(c())
  } else {
    return(c(parent,citrus.getClusterAncestors(parent,mergeOrder)))
  }
}

citrus.preCluster = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds=5,folds=NULL,transformCols=NULL,clusterConditionList=NULL,conditionComparaMatrix=NULL,balanceFactor=NULL,transformCofactor=5,...){
  
  addtlArgs = list(...)
  
  res=list()
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  if (!is.null(clusterConditionList)){
    allConditions = clusterConditionList
  } else if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
  } else {
    allConditions = as.list(colnames(fileList))
  }
  
  if (is.null(balanceFactor)){
    balanceFactor = as.factor(sample(c(0,1),nrow(fileList),replace=T))
  } else if (length(levels(balanceFactor))==1){
    balanceFactor = as.factor(sample(c(0,1),length(balanceFactor),replace=T))
  }
  if (nFolds=="all"){
    folds = list()
    nAllFolds=1
  } else if (!is.null(folds)){
    nAllFolds = length(folds)+1  
  } else {
    folds = pamr:::balanced.folds(y=balanceFactor,nfolds=nFolds)
    nAllFolds = nFolds+1  
  }
  folds[[nAllFolds]]="all"
  
  if (conditionCluster %in% names(addtlArgs)){
    clusterResult = parLapply(addtlArgs[["conditionCluster"]],allConditions,citrus.preClusterCondition,dataDir=dataDir,fileList=fileList,clusterCols=clusterCols,folds=folds,nFolds=nFolds,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,outputDir=outputDir,...)  
  } else {
    clusterResult = lapply(allConditions,citrus.preClusterCondition,dataDir=dataDir,fileList=fileList,clusterCols=clusterCols,folds=folds,nFolds=nFolds,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,outputDir=outputDir,...)  
  }
    
  names(clusterResult) = lapply(allConditions,paste,collapse="_vs_")
  
  save(clusterResult,file=file.path(outputDir,"citrus.clusterResult.rDat"))
  
  return(clusterResult)
}

citrus.preClusterCondition = function(conditions,dataDir,fileList,clusterCols,folds,nFolds,fileSampleSize,transformCols,transformCofactor,outputDir,...){
  
  cat(paste("Clustering Condition",paste(conditions,collapse=" vs "),"\n"))
  
  citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,...)
  
  foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions,...)
  cat("Assigning Events to Clusters\n")
  foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership,...)  
  
  preclusterObject = list(folds=folds,foldsCluster=foldsCluster,foldsClusterAssignments=foldsClusterAssignments,conditions=conditions,citrus.dataArray=citrus.dataArray,clusterColumns=clusterCols)
  if (nFolds!="all"){
    cat("Assigning Leftout Events to Clusters\n")
    leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols,...)  
    preclusterObject$leftoutClusterAssignments=leftoutClusterAssignments
  }
  
  return(preclusterObject)
  
}
  