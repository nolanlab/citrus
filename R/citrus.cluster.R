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

citrus.cluster = function(citrus.combinedFCSSet,clusteringColumns,clusteringType="hierarchical",...){
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

citrus.clusterFold = function(foldIndex,folds,citrus.combinedFCSSet,clusteringColumns,...){
  leavoutFileIds = as.vector(citrus.combinedFCSSet$fileIds[folds[[foldIndex]],])
  includeFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),leavoutFileIds)
  citrus.cluster(citrus.maskCombinedFCSSet(citrus.combinedFCSSet,fileIds=includeFileIds),clusteringColumns,...)
}

citrus.clusterAndMapFolds = function(citrus.combinedFCSSet,clusteringColumns,labels,nFolds=10,...){
  
  result = list()
  
  if (nFolds>1){
    # Define Folds
    result$folds = pamr:::balanced.folds(y=labels,nfolds=nFolds)
    
    # FOLDS MUST BE SORTED
    result$folds = lapply(result$folds,sort)
    
    # Cluster each fold of data together
    result$foldClustering = lapply(1:nFolds,citrus.clusterFold,folds=result$folds,citrus.combinedFCSSet=citrus.combinedFCSSet,clusteringColumns=clusteringColumns,...)
    
    # Map the left out data to the fold clustered data
    result$foldMappingAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,folds=result$folds,foldClustering=result$foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet)  
  }
  
  # Cluster all the data together
  result$allClustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns,...)
  
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


