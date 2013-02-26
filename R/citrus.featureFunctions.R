citrus.getFeatureSetNames = function(){
  return(c("densities","medians"))
}

citrus.buildFoldFeatures = function(index,featureTypes=c("densities"),folds,citrus.dataArray,foldsClusterAssignments,foldLargeEnoughClusters,conditions,calculateLeaveoutData=F,...){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions])
  } else if (calculateLeaveoutData){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
  } else {
    foldsFileIds=as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
  }
  return(citrus.buildFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],largeEnoughClusters=foldLargeEnoughClusters[[index]],foldsFileIds=foldsFileIds,foldFileNames=citrus.dataArray$fileNames[foldsFileIds],...))
}

citrus.buildFeatures = function(clusterAssignments,featureTypes,data,largeEnoughClusters,foldsFileIds,foldFileNames,...){
  addtlArgs = list(...)
  features = list()
  if ("densities" %in% featureTypes){
    cat(paste("Calculating Densities for files:",paste(foldsFileIds,collapse=", "),"\n"))
    densityFeatures = t(sapply(foldsFileIds,citrus.calculateFileClusterDensities,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,fileIds=data[,"fileId"]))
    rownames(densityFeatures) = foldFileNames
    colnames(densityFeatures) = paste("cluster",largeEnoughClusters,"density")
    features[["densityFeatures"]]=densityFeatures
  }
  if ("medians" %in% featureTypes){
    if (!("medianColumns" %in% names(addtlArgs))){
      stop("medianColumns argument must be specified to compute cluster medians.")
    }
    cat(paste("Calculating medians for files:",paste(foldsFileIds,collapse=", "),"\n"))
    medianFeatures = t(sapply(foldsFileIds,citrus.calculateFileClusterMedians,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,data=data,medianColumns=addtlArgs[["medianColumns"]]))
    rownames(medianFeatures) = foldFileNames
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
  
  # Are there zero cells assigned to our cluster?
  if (length(clusterAssignments[[clusterId]])==0){
    medians = rep(0,length(medianColumns))
  
  } else if (length(clusterAssignments[[clusterId]])==1) {
    # Are there 1 cells assigned to our cluster?
    
    if (fileId==include["fileId"]){
      medians = rep(0,length(medianColumns))
    } else {
      medians = include[medianColumns]
    }
    
  } else {
    # Are there more than 1 cells assigned to our cluster?
    
    include = include[include[,"fileId"]==fileId,medianColumns]
    if (length(medianColumns)>1){
      # do we have data?
      if (length(include)>0){
        if (is.null(nrow(include))){
          medians = include  
        } else {
          medians = apply(include,2,median) 
        }
      } else {
        medians = (rep(0,length(medianColumns)))
      }
    } else if (length(medianColumns)==1){
      if (length(include)>0){
        medians = median(include)
      } else {
        medians=0
      }
    }
  }
  
  names(medians) = paste(paste(paste("cluster",clusterId),medianColumns),"median")
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

citrus.getFeatureTypes = function(){
  return(c("densities","medians"))
}