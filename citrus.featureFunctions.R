citrus.calculateClusterFeatures = function(features="densities",
                                        minimumClusterSizePercent=0.05,
                                        citrus.dataArray,
                                        clusterConditionMatrix,
                                        folds,
                                        medianCols=NULL,
                                        workingDirectory){
  
  if (!file.exists(workingDirectory)){
    stop(paste("Directory",workingDirectory,"not found."))
  }
  conditions = rownames(clusterConditionMatrix)
  for (j in 1:ncol(clusterConditionMatrix)){
    for (i in 1:nrow(clusterConditionMatrix)){
      if (clusterConditionMatrix[i,j]){
        matrixIndex = (j-1)*length(conditions)+i
        load(paste(workingDirectory,"clusterChildren_",matrixIndex,".rDat",sep=""))
        load(file=paste(workingDirectory,"leaveoutChildren_",matrixIndex,".rDat",sep=""))
        
        if ("densities" %in% features){
          cat(paste("Calculating cluster densities densities for condition",paste(unique(conditions[c(i,j)]),collapse=" vs. "),"\n"))
          densityFeatures = lapply(as.list(1:length(clustered)),citrus.calculateClusterDensities,folds=folds,conditions=unique(c(conditions[i],conditions[j])),citrus.dataArray=citrus.dataArray,clusterChildren=clusterChildren,leaveoutChildren=leaveoutChildren,minimumClusterSizePercent=minimumClusterSizePercent)
          #citrus.calculateClusterDensities(1,folds=folds,conditions=unique(c(conditions[i],conditions[j])),citrus.dataArray=citrus.dataArray,clusterChildren=clusterChildren,leaveoutChildren=leaveoutChildren,minimumClusterSizePercent=minimumClusterSizePercent)
          if (conditions[i]!=conditions[j]){
            densityFeatures = lapply(densityFeatures,citrus.calculateConditionFeatureDifferences,condition1=citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[i]]],condition2=citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[j]]])
          }
          save(densityFeatures,file=paste(workingDirectory,"densityFeatures_",matrixIndex,".rDat",sep=""))
        }
        if ("clusterMedians" %in% features){
          cat(paste("Calculating cluster medians for condition",paste(unique(conditions[c(i,j)]),collapse=" vs. "),"\n"))
        }
        
      }
    }
  }
  
}

citrus.calculateConditionFeatureDifferences = function(data,condition1,condition2){
  x = data[condition2,]-data[condition1,]
  rownames(x)=apply(cbind(rownames(data[condition2,]),rownames(data[condition1,])),1,paste,collapse=" - ")
  return(x)
}

citrus.calculateClusterDensities = function(index,folds,conditions,citrus.dataArray,clusterChildren,leaveoutChildren,minimumClusterSizePercent){
  
  trainDensities = t(sapply(clusteredFileIds,citrus.calcFileClusterDensities,clusterIds=largeEnoughClusters,clusterChildren=clusterChildren[[index]],fileIds=clusteredDataFileIds)) 
  colnames(trainDensities) = paste("cluster",largeEnoughClusters,"density")
  rownames(trainDensities) = citrus.dataArray$fileNames[clusteredFileIds]
  
  if (!(index > length(folds))){
    testDensities = t(sapply(leaveoutFileIds,citrus.calcFileClusterDensities,clusterIds=largeEnoughClusters,clusterChildren=leaveoutChildren[[index]],fileIds=leaveoutDataFileIds))
    colnames(testDensities) = paste("cluster",largeEnoughClusters,"density")
    rownames(testDensities) = citrus.dataArray$fileNames[leaveoutFileIds]
    combined = rbind(trainDensities,testDensities)
    return(combined[sapply(citrus.dataArray$fileNames[,conditions],citrus.indexof,y=rownames(combined)),])
  } else {
    return(trainDensities)
  }
}

citrus.calcFileClusterDensities = function(fileId,clusterIds,clusterChildren,fileIds){
  sapply(clusterIds,citrus.calcFileClusterDensity,clusterChildren=clusterChildren,fileId=fileId,fileIds=fileIds)
}

citrus.calcFileClusterDensity = function(clusterId,clusterChildren,fileId,fileIds){
  sum(which(fileIds==fileId) %in% clusterChildren[[clusterId]])/sum((fileIds==fileId))
}

citrus.buildFoldFeatures = function(index,featureTypes=c("densities"),citrus.dataArray,foldsClusterAssignments,foldLargeEnoughClusters,conditions,minimumClusterSizePercent=0.05,calculateLeaveoutData=F){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions])
  } else if (calculateLeaveoutData){
    foldsFileIds=as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
  } else {
    foldsFileIds=as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
  }
  return(citrus.buildFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,citrus.dataArray=citrus.dataArray,largeEnoughClusters=foldLargeEnoughClusters[[index]],fileIds=foldsFileIds))
}

citrus.buildFeatures = function(clusterAssignments,featureTypes,citrus.dataArray,largeEnoughClusters,fileIds){
  if ("densities" %in% featureTypes){
    densityFeatures = t(sapply(fileIds,citrus.calculateFileClusterDensities,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,fileIds=citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%fileIds,"fileId"]))
    rownames(densityFeatures) = citrus.dataArray$fileNames[fileIds]
    colnames(densityFeatures) = paste("cluster",largeEnoughClusters,"density")
    return(densityFeatures)
  }
}

citrus.calculateFileClusterDensities = function(fileId,clusterIds,clusterAssignments,fileIds){
  sapply(clusterIds,citrus.calculateFileClusterDensity,clusterAssignments=clusterAssignments,fileId=fileId,fileIds=fileIds)
}

citrus.calculateFileClusterDensity = function(clusterId,clusterAssignments,fileId,fileIds){
    sum(which(fileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((fileIds==fileId))
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


