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
  return(citrus.buildFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,largeEnoughClusters=foldLargeEnoughClusters[[index]],foldsFileIds=foldsFileIds,conditions=conditions,citrus.dataArray=citrus.dataArray,...))
}

citrus.buildFeatures = function(clusterAssignments,featureTypes,largeEnoughClusters,foldsFileIds,conditions,citrus.dataArray,...){
  features = list()
  for (featureType in featureTypes){
    #features[[featureType]]=do.call(paste("citrus.calculateFeature",featureType,sep="."),args=list(foldsFileIds=foldsFileIds,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],conditions=conditions,citrus.dataArray=citrus.dataArray,emdColumns=emdColumns))
    features[[featureType]] = do.call(paste("citrus.calculateFeature",featureType,sep="."),args=list(foldsFileIds=foldsFileIds,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],conditions=conditions,citrus.dataArray=citrus.dataArray,...))
    
    # ASSUME THAT WANT FEATURE DIFFERENCES. MAY BE BAD ASSUMPTION
    if ((length(conditions)==2) && (featureType!="emDists")){
      cat("Two conditions found. Assuming differential features of interest.\n")
      cat(paste("Caluclating difference in ",featureType," between ",conditions[2]," & ",conditions[1],".\n",sep=""))
      fns2 = citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[2]]]
      fns1 = citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[1]]]
      features[[featureType]] = features[[featureType]][rownames(features[[featureType]]) %in% fns2,] - features[[featureType]][rownames(features[[featureType]]) %in% fns1,]
      colnames(features[[featureType]]) = paste(colnames(features[[featureType]]),"difference")
    }
    
  }
  return(do.call("cbind",features))
}

citrus.calculateFeature.emDists = function(foldsFileIds,clusterIds,clusterAssignments,data,foldFileNames,conditions,citrus.dataArray,...){
  addtlArgs = list(...)
  if (!("emdColumns" %in% names(addtlArgs))){
    stop("emdColumns argument must be specified to compute cluster emDists.")
  }
  if (length(conditions)!=2){
    stop("Only know how to calculate EMD Features for two conditions.")
  }
  referenceFileIds = foldsFileIds[foldsFileIds %in% citrus.dataArray$fileIds[,conditions[1]]]
  targetFileIds = foldsFileIds[foldsFileIds %in% citrus.dataArray$fileIds[,conditions[2]]]
  #features = t(sapply(1:length(referenceFileIds),citrus.calculateFileClustersEMDist,clusterIds=clusterIds,clusterAssignments=clusterAssignments,referenceFileIds=referenceFileIds,targetFileIds=targetFileIds,data=data,emdColumns=emdColumns))
  features = t(sapply(1:length(referenceFileIds),citrus.calculateFileClustersEMDist,clusterIds=clusterIds,clusterAssignments=clusterAssignments,referenceFileIds=referenceFileIds,targetFileIds=targetFileIds,data=data,emdColumns=addtlArgs[["emdColumns"]]))
  rownames(features) = citrus.dataArray$fileNames[targetFileIds]
  return(features)
}

citrus.calculateFileClustersEMDist = function(sampleIndex,clusterIds,clusterAssignments,referenceFileIds,targetFileIds,data,emdColumns){
  res = sapply(clusterIds,citrus.calculateFileClusterEMDist,clusterAssignments=clusterAssignments,referenceFileId=referenceFileIds[sampleIndex],targetFileId=targetFileIds[sampleIndex],data=data,emdColumns=emdColumns)
  res2 = as.vector(res) 
  names(res2) = paste(paste(paste("cluster",rep(clusterIds,each=nrow(res))),rownames(res)),"emDist")
  return(res2)  
}

citrus.calculateFileClusterEMDist = function(clusterId,clusterAssignments,referenceFileId,targetFileId,data,emdColumns){
  clusterData = data[clusterAssignments[[clusterId]],]
  referenceData = clusterData[clusterData[,"fileId"]==referenceFileId,]
  targetData = clusterData[clusterData[,"fileId"]==targetFileId,]
  sapply(emdColumns,citrus.calculateFileClusterParameterEMDist,referenceData,targetData)
}

citrus.calculateFileClusterParameterEMDist = function(emdColumn,referenceData,targetData){
  cat("IMPLEMENT MINIMUM CLUSTER PERCENTAGE CHECK\n");
  h = hist(c(referenceData[,emdColumn],targetData[,emdColumn]),breaks=50,plot=F)
  emdw(A=h$mids,wA=hist(referenceData[,emdColumn],plot=F,breaks=h$breaks)$density,B=h$mids,wB=hist(targetData[,emdColumn],plot=F,breaks=h$breaks)$density)
}

citrus.calculateFeature.densities = function(foldsFileIds,clusterIds,clusterAssignments,data,citrus.dataArray,...){
  features = t(sapply(foldsFileIds,citrus.calculateFileClustersDensities,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))
  rownames(features) = citrus.dataArray$fileNames[foldsFileIds]
  return(features)
}

citrus.calculateFileClustersDensities = function(fileId,clusterIds,clusterAssignments,data,...){
  fileIds=data[,"fileId"]
  res = sapply(clusterIds,citrus.calculateFileClusterDensity,clusterAssignments=clusterAssignments,fileId=fileId,fileIds=fileIds)
  names(res) = paste(paste("cluster",clusterIds),"density")
  return(res)
}

citrus.calculateFileClusterDensity = function(clusterId,clusterAssignments,fileId,fileIds){
  sum(which(fileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((fileIds==fileId))
}

citrus.calculateFeature.medians = function(foldsFileIds,clusterIds,clusterAssignments,data,citrus.dataArray,...){
  features = t(sapply(foldsFileIds,citrus.calculateFileClustersMedians,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))
  rownames(features) = citrus.dataArray$fileNames[foldsFileIds]
  return(features)
}

citrus.calculateFileClustersMedians = function(fileId,clusterIds,clusterAssignments,data,...){
  addtlArgs = list(...)
  if (!("medianColumns" %in% names(addtlArgs))){
    stop("medianColumns argument must be specified to compute cluster medians.")
  }
  unlist(lapply(clusterIds,citrus.calculateFileClusterMedian,clusterAssignments=clusterAssignments,fileId=fileId,data=data,medianColumns=addtlArgs[["medianColumns"]]))
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

