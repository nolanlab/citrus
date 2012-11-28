citrus.buildFoldFeatures = function(index,featureTypes=c("densities"),citrus.dataArray,foldsClusterAssignments,foldLargeEnoughClusters,conditions,calculateLeaveoutData=F,...){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions])
  } else if (calculateLeaveoutData){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
  } else {
    foldsFileIds=as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
  }
  return(citrus.buildFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],largeEnoughClusters=foldLargeEnoughClusters[[index]],...))
}

#citrus.buildFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],largeEnoughClusters=foldLargeEnoughClusters[[index]])

citrus.buildFeatures = function(clusterAssignments,featureTypes,data,largeEnoughClusters,...){
  addtlArgs = list(...)
  features = list()
  foldFileIds = unique(data[,"fileId"])
  if ("densities" %in% featureTypes){
    cat(paste("Calculating Densities for files:",paste(foldFileIds,collapse=", "),"\n"))
    densityFeatures = t(sapply(foldFileIds,citrus.calculateFileClusterDensities,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,fileIds=data[,"fileId"]))
    rownames(densityFeatures) = citrus.dataArray$fileNames[foldFileIds]
    colnames(densityFeatures) = paste("cluster",largeEnoughClusters,"density")
    singularClusterDensity = which(apply(densityFeatures==1,2,all))
    if (length(singularClusterDensity)>0){
      densityFeatures = densityFeatures[,-singularClusterDensity]
    }
    features[["densityFeatures"]]=densityFeatures
  }
  if ("medians" %in% featureTypes){
    if (!("medianColumns" %in% names(addtlArgs))){
      stop("medianColumns argument must be specified to compute cluster medians.")
    }
    cat(paste("Calculating medians for files:",paste(foldFileIds,collapse=", "),"\n"))
    medianFeatures = t(sapply(foldFileIds,citrus.calculateFileClusterMedians,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,data=data,medianColumns=addtlArgs[["medianColumns"]]))
    rownames(medianFeatures) = citrus.dataArray$fileNames[foldFileIds]
    features[["medianFeatures"]]=medianFeatures
  }
  return(do.call("cbind",features))
}

citrus.calculateFileClusterDensities = function(fileId,clusterIds,clusterAssignments,fileIds){
  sapply(clusterIds,citrus.calculateFileClusterDensity,clusterAssignments=clusterAssignments,fileId=fileId,fileIds=fileIds)
}

citrus.calculateFileClusterDensity = function(clusterId,clusterAssignments,fileId,fileIds){
  sum(which(fileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((fileIds==fileId))
}

citrus.calculateFileClusterMedians = function(fileId,clusterIds,clusterAssignments,data,medianColumns){
  unlist(lapply(clusterIds,citrus.calculateFileClusterMedian,clusterAssignments=clusterAssignments,fileId=fileId,data=data,medianColumns=medianColumns))
}

citrus.calculateFileClusterMedian = function(clusterId,clusterAssignments,fileId,data,medianColumns){
  include = data[clusterAssignments[[clusterId]],]
  include = include[include[,"fileId"]==fileId,medianColumns]
  if (is.null(nrow(include))){
    if (length(include)>0){
      return(include)
    } else {
      # Assume empty clusters have 0 values. Maybe a bad assumption. 
      return(rep(0,length(medianColumns)))
    }
  }
  medians = apply(include,2,median)
  names(medians) = paste(paste(paste("cluster",clusterId),colnames(data)[medianColumns]),"median")
  return(medians)
}


citrus.calculateConditionFeatureDifferences = function(data,condition1,condition2){
  x = data[condition2,]-data[condition1,]
  rownames(x)=apply(cbind(rownames(data[condition2,]),rownames(data[condition1,])),1,paste,collapse=" - ")
  return(x)
}

citrus.calculateFoldLargeEnoughClusters = function(index,foldsClusterAssignments,folds,citrus.dataArray,conditions,minimumClusterSizePercent=0.05){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    foldsDataLength = sum(citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions])
  } else {
    foldsDataLength = sum(citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[-folds[[index]],conditions])
  }
  minimumClusterSize = floor(minimumClusterSizePercent*foldsDataLength)
  return(citrus.calculateLargeEnoughClusters(foldsClusterAssignments[[index]],minimumClusterSize))
}
  
citrus.calculateLargeEnoughClusters = function(clusters,minimumClusterSize){
  clusterLengths = lapply(clusters,length)
  return(which(clusterLengths >= minimumClusterSize))
}
