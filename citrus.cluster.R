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
  nnMap = SPADE.assignToCluster(tbl=newData,cluster_data=data,cluster_assign=rep(1,nrow(data)))
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
