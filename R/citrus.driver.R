citrus.full = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,transformCols=NULL,...){
  
  # Error check before we actually start the work.
  if ((!all(featureTypes %in% citrus.getFeatureTypes()))||(length(featureTypes)<1)){
    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.getFeatureTypes(),collapse=", "),"\n"))
  }
  
  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  
  labelCol = which(colnames(fileList)=="labels")
  if (length(labelCol)==0){
    stop("'labels' column not found in file list.")
  }
  
  if ("conditionComparaMatrix" %in% names(list(...))){
    #allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
    allConditions = citrus.convertConditionMatrix(list(...)[["conditionComparaMatrix"]]) 
  } else {
    allConditions = as.list(colnames(fileList)[-labelCol])
  }

  for (conditions in allConditions){
    cat(paste("Analyzing Condition",paste(conditions,collapse=" vs "),"\n"))
    
    conditionOutputDir = file.path(outputDir,paste(conditions,collapse="_vs_"))
    
    citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize)
    
    
    nAllFolds = nFolds+1
    folds = pamr:::balanced.folds(y=fileList[,"labels"],nfolds=nAllFolds)
    folds[[nAllFolds]]="all"
      
    
    cat("Clustering\n")
    foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions)
    cat("Assigning Events to Clusters\n")
    foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)
    cat("Assigning Leftout Events to Clusters\n")
    leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols)
    cat("Calculating Fold Large Enough Clusters\n")
    foldLargeEnoughClusters = lapply(1:nAllFolds,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=foldsClusterAssignments,folds=folds,citrus.dataArray=citrus.dataArray,minimumClusterSizePercent=0.05)
      
    foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,...)
    #foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=medianColumns)
    leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,...)
    #leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,medianColumns=medianColumns)
    
    if (length(conditions)==2){
      foldFeatures = lapply(1:length(folds),citrus.buildFoldFeatureDifferences,features=foldFeatures,conditions=conditions,citrus.dataArray=citrus.dataArray,folds=folds)
      leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatureDifferences,features=leftoutFeatures,conditions=conditions,citrus.dataArray=citrus.dataArray,folds=folds,calculateLeaveoutData=T)
    }
    
    regularizationThresholds = citrus.generateRgularizationThresholds(foldFeatures[[nAllFolds]],labels=fileList[,"labels"],modelTypes=modelTypes)
    
    foldModels = lapply(modelTypes,citrus.buildTypeModels,folds=folds,foldFeatures=foldFeatures,labels=fileList[,"labels"],regularizationThresholds=regularizationThresholds)
    names(foldModels)=modelTypes
    
    leftoutPredictions = lapply(modelTypes,citrus.foldTypePredict,foldModels=foldModels,leftoutFeatures=leftoutFeatures)
    names(leftoutPredictions)=modelTypes
    
    predictionSuccess = lapply(as.list(modelTypes),citrus.foldTypeScore,folds=folds,leftoutPredictions=leftoutPredictions,labels=fileList[,"labels"])
    names(predictionSuccess)=modelTypes
    
    thresholdSEMs = lapply(modelTypes,citrus.modelTypeSEM,predictionSuccess=predictionSuccess)
    names(thresholdSEMs)=modelTypes
    
    thresholdErrorRates = lapply(modelTypes,citrus.calcualteTypeErroRate,predictionSuccess=predictionSuccess)
    names(thresholdErrorRates)=modelTypes
    
    thresholdFDRRates = lapply(modelTypes,citrus.calculateTypeFDRRate,foldModels=foldModels,foldFeatures=foldFeatures,labels=fileList[,"labels"])
    names(thresholdFDRRates)=modelTypes
    
    cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdErrorRates=thresholdErrorRates,thresholdSEMs=thresholdSEMs,thresholdFDRRates=thresholdFDRRates)
    names(cvMinima)=modelTypes
    
    # Make condition output directoy
    dir.create(conditionOutputDir,showWarnings=T,recursive=T)
    
    # Plot
    sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=regularizationThresholds,thresholdErrorRates=thresholdErrorRates,thresholdFDRRates=thresholdFDRRates,cvMinima=cvMinima,thresholdSEMs=thresholdSEMs,foldModels=foldModels)
    
    # Extract Features
    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
    #citrus.extractModelFeatures("pamr",cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
    names(differentialFeatures) = modelTypes
    
    # Plot Features
    lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=differentialFeatures,foldFeatures=foldFeatures,outputDir=conditionOutputDir,labels=fileList[,"labels"])
    
    # Plot Clusters
    lapply(modelTypes,citrus.plotClusters,differentialFeatures=differentialFeatures,outputDir=conditionOutputDir,clusterChildren=foldsClusterAssignments,citrus.dataArray=citrus.dataArray,conditions=conditions,clusterCols=clusterCols)
  }
}


citrus.buildFoldFeatureDifferences = function(index,features,conditions,citrus.dataArray,folds,calculateLeaveoutData=F){
  
  if (folds[[index]]=="all"){
    includeRowIds=1:nrow(citrus.dataArray$fileIds)
  } else if (!calculateLeaveoutData){
    includeRowIds=setdiff(1:nrow(citrus.dataArray$fileIds),folds[[index]])
  } else {
    includeRowIds = folds[[index]]
  }
  diffFeatures = features[[index]][citrus.dataArray$fileNames[citrus.dataArray$fileIds[includeRowIds,conditions[2]]],] - features[[index]][citrus.dataArray$fileNames[citrus.dataArray$fileIds[includeRowIds,conditions[1]]],] 
  colnames(diffFeatures) = paste(colnames(diffFeatures),"difference")
  return(diffFeatures)
}
