#####################
# Mapping new Stuff
#####################

citrus.mapToClusterSpace = function(citrus.combinedFCSSet.new,citrus.combinedFCSSet.old,citrus.clustering,mappingColumns=NULL,...){
  if (is.null(mappingColumns)){
    mappingColumns = citrus.clustering$clusteringColumns
  }
  
  # Map new data to nearest neighbor in existing clustered dataset
  cat(paste("Mapping",nrow(citrus.combinedFCSSet.new$data),"events\n"))
  nearestNeighborMap = citrus.assignToCluster(tbl=citrus.combinedFCSSet.new$data[,mappingColumns],cluster_data=citrus.combinedFCSSet.old$data[,mappingColumns],cluster_assign=rep(1,nrow(citrus.combinedFCSSet.old$data)))
  
  # Assign new data to the same clusters as its nearest neighbor in the existing clustered dataset
  newDataClusterAssignments = mclapply(citrus.clustering$clusterMembership,citrus.mapNeighborsToCluster,nearestNeighborMap=nearestNeighborMap,...)
  
  result = list(clusterMembership=newDataClusterAssignments,mappingColumns=mappingColumns,call=match.call())
  class(result) = "citrus.mapping"
  return(result)
}

citrus.mapNeighborsToCluster = function(clusterMembers,nearestNeighborMap){
  which(nearestNeighborMap %in% clusterMembers)
}

citrus.mapFoldDataToClusterSpace = function(foldIndex,folds,foldClustering,citrus.combinedFCSSet,...){
  leavoutFileIds = as.vector(citrus.combinedFCSSet$fileIds[folds[[foldIndex]],])
  includeFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),leavoutFileIds)
  
  
  citrus.mapToClusterSpace(citrus.combinedFCSSet.new=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,leavoutFileIds),
                           citrus.combinedFCSSet.old=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,includeFileIds),
                           citrus.clustering=foldClustering[[foldIndex]],
                           mappingColumns=foldClustering[[foldIndex]]$clusteringColumns,...)
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

citrus.cluster = function(citrus.combinedFCSSet,clusteringColumns,clusteringType="hierarchical"){
  clustering = do.call(paste0("citrus.cluster.",clusteringType),args=list(data=citrus.combinedFCSSet$data[,clusteringColumns]))
  clusterMembership = do.call(paste0("citrus.calculateClusteringMembership.",clusteringType),args=list(clustering=clustering))
  result = list(clustering=clustering,clusterMembership=clusterMembership,type=clusteringType,clusteringColumns=clusteringColumns,call=match.call())
  class(result) = "citrus.clustering"
  return(result)
}

print.citrus.clustering = function(x,...){
  cat("Citrus clustering\n")
  cat(paste0("\tType:\t\t",x$type,"\n"))
  cat(paste0("\tClusters:\t",length(x$clusterMembership),"\n"))
}

citrus.clusterFold = function(foldIndex,folds,citrus.combinedFCSSet,clusteringColumns,clusteringType="hierarchical"){
  leavoutFileIds = as.vector(citrus.combinedFCSSet$fileIds[folds[[foldIndex]],])
  includeFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),leavoutFileIds)
  citrus.cluster(citrus.maskCombinedFCSSet(citrus.combinedFCSSet,fileIds=includeFileIds),clusteringColumns,clusteringType)
}

citrus.clusterAndMapFolds = function(citrus.combinedFCSSet,clusteringColumns,labels,clusteringType="hierarchical",nFolds=10,...){
  
  result = list()
  
  if (nFolds>1){
    # Define Folds
    result$folds = pamr:::balanced.folds(y=labels,nfolds=nFolds)  
    
    # Cluster each fold of data together
    result$foldClustering = lapply(1:nFolds,citrus.clusterFold,folds=result$folds,citrus.combinedFCSSet=citrus.combinedFCSSet,clusteringColumns=clusteringColumns,clusteringType=clusteringType)
    
    # Map the left out data to the fold clustered data
    result$foldMappingAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,folds=result$folds,foldClustering=result$foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet)  
  }
  
  # Cluster all the data together
  result$allClustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns,clusteringType)
  
  result$nFolds = nFolds
  
  class(result) = "citrus.foldClustering"
  return(result)
}

print.citrus.foldClustering = function(x,...){
  cat("citrus.foldClustering\n")
  cat(paste("\tNumber of Folds: ",length(x$folds),"\n"))
}



citrus.cluster.hierarchical = function(data){
  cat(paste("Clustering",nrow(data),"events\n"));
  return(Rclusterpp.hclust(data))  
}

citrus.calculateClusteringMembership.hierarchical = function(clustering,...){
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

#citrus.preCluster = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds=5,folds=NULL,transformCols=NULL,clusterConditionList=NULL,conditionComparaMatrix=NULL,balanceFactor=NULL,transformCofactor=5,...){
#  
#  addtlArgs = list(...)
#  
#  res=list()
#  if (!file.exists(outputDir)){
#    stop(paste("Output directory",outputDir,"not found."))
#  }
#  
#  if (!is.null(clusterConditionList)){
#    allConditions = clusterConditionList
#  } else if (!is.null(conditionComparaMatrix)){
#    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
#  } else {
#    allConditions = as.list(colnames(fileList))
#  }
#  
#  if (is.null(balanceFactor)){
#    balanceFactor = as.factor(sample(c(0,1),nrow(fileList),replace=T))
#  } else if (length(levels(balanceFactor))==1){
#    balanceFactor = as.factor(sample(c(0,1),length(balanceFactor),replace=T))
#  }
#  if (nFolds=="all"){
#    folds = list()
#    nAllFolds=1
#  } else if (!is.null(folds)){
#    nAllFolds = length(folds)+1  
#  } else {
#    folds = pamr:::balanced.folds(y=balanceFactor,nfolds=nFolds)
#    nAllFolds = nFolds+1  
#  }
#  folds[[nAllFolds]]="all"
#  
#  
#  if ("conditionParallelClusters" %in% names(addtlArgs)){
#    clusterResult = parLapply(addtlArgs[["conditionParallelClusters"]],allConditions,citrus.preClusterCondition,dataDir=dataDir,fileList=fileList,clusterCols=clusterCols,folds=folds,nFolds=nFolds,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,outputDir=outputDir,...)  
#  } else {
#    clusterResult = lapply(allConditions,citrus.preClusterCondition,dataDir=dataDir,fileList=fileList,clusterCols=clusterCols,folds=folds,nFolds=nFolds,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,outputDir=outputDir,...)  
#  }
#    
#  names(clusterResult) = lapply(allConditions,paste,collapse="_vs_")
#  
#  save(clusterResult,file=file.path(outputDir,"citrus.clusterResult.rDat"))
#  
#  return(clusterResult)
#}

#citrus.preClusterCondition = function(conditions,dataDir,fileList,clusterCols,folds,nFolds,fileSampleSize,transformCols,transformCofactor,outputDir,...){
  
#  cat(paste("Clustering Condition",paste(conditions,collapse=" vs "),"\n"))
  
#  addtlArgs = list(...)
  
#  citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,fileSampleSize=fileSampleSize,transformCols=transformCols,transformCofactor=transformCofactor,...)
  
#  foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions,...)
#  cat("Assigning Events to Clusters\n")
#  foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership,...)  
#  
#  preclusterObject = list(folds=folds,foldsCluster=foldsCluster,foldsClusterAssignments=foldsClusterAssignments,conditions=conditions,citrus.dataArray=citrus.dataArray,clusterColumns=clusterCols)
#  if (nFolds!="all"){
#    cat("Assigning Leftout Events to Clusters\n")
#    leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols,...)  
#    preclusterObject$leftoutClusterAssignments=leftoutClusterAssignments
#  }
#  
#  return(preclusterObject)
  
#}
  

#


#citrus.mapDataToClusterSpace = function(dataDir,newFileList,preClusterResult,fileSampleSize,mappingColumns=NULL,transformCols=NULL,transformCofactor=5,...){
#  conditions = strsplit(names(preClusterResult),"_vs_")
#  mappingResult = list()
#  for (conditionName in names(preClusterResult)){
#    conditions = strsplit(conditionName,"_vs_")[[1]]
#    mappingResult[[conditionName]]$citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=newFileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize,transformCofactor=transformCofactor,...)
#    
#    if (is.null(mappingColumns)){
#      mappingColumns = preClusterResult[[conditionName]]$clusterColumns
#    }
#    
#    cat(paste("Mapping",nrow(mappingResult[[conditionName]]$citrus.dataArray$data),"events to existing cluster space.\n"));
#    mappingResult[[conditionName]]$foldsClusterAssignments = list()
#    mappingResult[[conditionName]]$foldsClusterAssignments[["mappingResult"]] = citrus.mapDataToClusterSpace(data=preClusterResult[[conditionName]]$citrus.dataArray$data[,mappingColumns],
#                                 clusterAssignments=preClusterResult[[conditionName]]$foldsClusterAssignments[[which(preClusterResult[[conditionName]]$folds=="all")]],
#                                 newData=mappingResult[[conditionName]]$citrus.dataArray$data[,mappingColumns],...)
#    mappingResult[[conditionName]]$folds = list("mappingResult")
#    mappingResult[[conditionName]]$conditions = conditions
#    mappingResult[[conditionName]]$clusterColumns = preClusterResult[[conditionName]]$clusterColumns
#  }
#  return(mappingResult)
#}


#citrus.mapFoldDataToClusterSpace = function(index,citrus.dataArray,foldClusterAssignments,folds,conditions,clusterCols,...){
#  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
#    return(NULL)
#  }
#  foldFileIds = as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
#  leftoutFileIds = as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
#  return(citrus.mapDataToClusterSpace(data=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%foldFileIds,clusterCols],clusterAssignments=foldClusterAssignments[[index]],newData=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%leftoutFileIds,clusterCols],...))  
#}

#citrus.mapDataToClusterSpace = function(data,clusterAssignments,newData,...){  
#  nnMap = citrus.assignToCluster(tbl=newData,cluster_data=data,cluster_assign=rep(1,nrow(data)))
#  addtlArgs = list(...)
#  if ("mc.cores" %in% names(addtlArgs)){
#    return(mclapply(clusterAssignments,citrus.mapNeighborsToCluster,nearestNeighborMap=nnMap,mc.cores=addtlArgs[["mc.cores"]]))
#  } else {
#    return(lapply(clusterAssignments,citrus.mapNeighborsToCluster,nearestNeighborMap=nnMap))
#  }
#}  

#citrus.mapNeighborsToCluster = function(clusterMembers,nearestNeighborMap){
#  which(nearestNeighborMap %in% clusterMembers)
#}

# Deal with 

#citrus.foldCluster = function(foldMembership,citrus.dataArray,clusterCols,conditions,...){
#  if ((length(foldMembership)==1) && (foldMembership=="all")){
#    cat("Clustering All\n")
#    includeFileIds = as.vector(citrus.dataArray$fileIds[,conditions])
#  } else {
#    includeFileIds = as.vector(citrus.dataArray$fileIds[-foldMembership,conditions])
#    cat(paste("Clustering fileIds",paste(includeFileIds,collapse=", "),"\n"))
#  }
#  if (any(!is.numeric(clusterCols))){
#    containedCols = setdiff(clusterCols,colnames(citrus.dataArray$data))
#    if (length(containedCols)>0){
#      stop(paste("Cluster cols",paste(containedCols,collapse=", "),"not found. Valid channel names:",paste(colnames(citrus.dataArray$data),collapse=", ")))
#    }  
#  }
#  return(citrus.cluster(citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% includeFileIds,clusterCols],...))
#}

#citrus.cluster = function(data,...){
#  cat(paste("Clustering",nrow(data),"events\n"));
#  
#  clusterType="hierarchical"
#  if ("clusterType" %in% names(list(...))){
#    clusterType = list(...)[["clusterType"]]
#  }
#  if (clusterType=="hierarchical"){
#    return(Rclusterpp.hclust(data))  
#  } else {
#    stop(paste("Don't know how to cluster using",clusterType))
#  }
#}


