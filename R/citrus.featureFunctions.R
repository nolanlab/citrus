###############################################
# Generic function for building features
###############################################
citrus.buildFeatures = function(citrus.combinedFCSSet,clusterAssignments,clusterIds,featureType="abundances",...){  

  #Error check before we actually start the work.
  if ((!all(featureType %in% citrus.featureTypes()))||(length(featureType)<1)){
    stop(paste("featureType must be one of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
  }
  
  do.call(paste0("citrus.calculateFeature.",featureType),args=list(clusterIds=clusterIds,clusterAssignments=clusterAssignments,citrus.combinedFCSSet=citrus.combinedFCSSet,...))
}

################################
# Abundance Features
################################
citrus.calculateFeature.abundances = function(clusterIds,clusterAssignments,citrus.combinedFCSSet,...){
  eventFileIds = citrus.combinedFCSSet$data[,"fileId"]
  fileIds = unique(eventFileIds)
  res = mcmapply(citrus.calculateFileClusterAbundance,
          clusterId=rep(clusterIds,length(fileIds)),
          fileId=rep(fileIds,each=length(clusterIds)),
          MoreArgs=list(
                        clusterAssignments=clusterAssignments,
                        eventFileIds=eventFileIds)
          ,...)
  return(
    matrix(res,ncol=length(clusterIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames,paste("cluster",clusterIds,"abundance")))
  )
}

citrus.calculateFileClusterAbundance = function(clusterId,fileId,clusterAssignments,eventFileIds){
  sum(which(eventFileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((eventFileIds==fileId))
}

################################
# Median Features
################################
citrus.calculateFeature.medians = function(clusterIds,clusterAssignments,citrus.combinedFCSSet,...){
  if (!("medianColumns" %in% names(list(...)))){
    stop("medianColumns argument missing")
  }
  
  medianColumns = list(...)[["medianColumns"]]
  fileIds = unique(citrus.combinedFCSSet$data[,"fileId"])
  res = mcmapply(citrus.calculateFileClusterMedian,
                 clusterId=rep(rep(clusterIds,length(medianColumns)),length(fileIds)),
                 fileId=rep(fileIds,each=(length(clusterIds)*length(medianColumns))),
                 medianColumn=rep(rep(medianColumns,each=length(clusterIds)),length(fileIds)),
                 MoreArgs=list(
                   data=citrus.combinedFCSSet$data,
                   clusterAssignments=clusterAssignments))
  
  return(
    matrix(res,nrow=length(fileIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames,paste("cluster",rep(clusterIds,length(medianColumns)),rep(medianColumns,each=length(clusterIds)),"median")))
  )
}

citrus.calculateFileClusterMedian = function(clusterId,fileId,medianColumn,data,clusterAssignments){
  clusterData = data[clusterAssignments[[clusterId]],]
  clusterFileDataValues = clusterData[clusterData[,"fileId"]==fileId,medianColumn]
  if (length(clusterFileDataValues)<3){
    return(0)
  } else {
    return(median(clusterFileDataValues))
  }
}


###############################################
# Functions for dealing with fold features
###############################################
citrus.buildFoldFeatureSet = function(citrus.foldClustering,citrus.combinedFCSSet,featureType="abundances",minimumClusterSizePercent=0.05,...){
  result = list()
  if (citrus.foldClustering$nFolds>1){
    # Select clusters in each folds
    # Default is minimum cluster size
    result$foldLargeEnoughClusters = lapply(citrus.foldClustering$foldClustering,citrus.selectClusters,minimumClusterSizePercent=minimumClusterSizePercent)  
    
    # Build Training Features
    result$foldFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.buildFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,foldClustering=citrus.foldClustering$foldClustering,...)
    
    # Build Testing Features
    result$leftoutFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.buildFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,foldMappingAssignments=citrus.foldClustering$foldMappingAssignments,calculateLeftoutFeatureValues=T,...)
  }
  
  # Build features for clustering of all samples
  result$allLargeEnoughClusters = citrus.selectClusters(citrus.clustering=citrus.foldClustering$allClustering,minimumClusterSizePercent=minimumClusterSizePercent)
  result$allFeatures = citrus.buildFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.foldClustering$allClustering$clusterMembership,clusterIds=result$allLargeEnoughClusters,featureType=featureType,...)
  
  # Extra feature building parameters, etc
  result$minimumClusterSizePercent=minimumClusterSizePercent
  result$folds=citrus.foldClustering$folds
  result$nFolds=citrus.foldClustering$nFolds
  result$minimumClusterSizePercent=minimumClusterSizePercent
  
  class(result) = "citrus.foldFeatureSet"
  return(result)
}

citrus.buildFoldFeatures = function(foldIndex,folds,foldClusterIds,citrus.combinedFCSSet,foldClustering=NULL,foldMappingAssignments=NULL,featureType="abundances",calculateLeftoutFeatureValues=F,...){
  leavoutFileIndices = folds[[foldIndex]]
  
  if (calculateLeftoutFeatureValues){
    featureFileIds = as.vector(citrus.combinedFCSSet$fileIds[leavoutFileIndices,])
    if (is.null(foldMappingAssignments)){
      stop("foldMappingAssignments argument missing")
    }
    cellAssignments = foldMappingAssignments[[foldIndex]]$clusterMembership
  } else {
    featureFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),as.vector(citrus.combinedFCSSet$fileIds[leavoutFileIndices,]))
    if (is.null(foldClustering)){
      stop("foldClustering argument missing")
    }
    cellAssignments = foldClustering[[foldIndex]]$clusterMembership
  }
  
  citrus.buildFeatures(citrus.combinedFCSSet=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,featureFileIds),
                       clusterAssignments=cellAssignments,
                       clusterIds=foldClusterIds[[foldIndex]],
                       featureType=featureType,
                       ...)
  
}