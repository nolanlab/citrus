citrus.full = function(dataDirectory,
                       clusteringColumns,
                       labels,
                       family,
                       fileList,
                       outputDirectory,
                       modelTypes=c("glmnet"),
                       nFolds=1,
                       plot=T,
                       conditions=NULL,
                       conditionComparaMatrix=NULL,
                       ...){
  
  # No point in running cv if SAM only model
  if (all(modelTypes=="sam")){
    nFolds=1
  }
  
  if (!is.null(conditions)){
    allConditions = as.list(conditions)
  } else if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix)
  } else {
    allConditions = as.list(colnames(fileList))
  }
  
  # Read in data
  citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory=dataDirectory,fileList=fileList,useChannelDescriptions=T,...)
    
  # Cluster each fold
  citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels=labels,nFolds=nFolds)
  save(list=c("citrus.combinedFCSSet","citrus.foldClustering"),file=file.path(outputDirectory,"citrusClustering.rData"),compress=F)
    
  conditionFeatures = list()
  conditionRegressionResults = list()
  # Analyze each condition or set of conditions separately 
  for (conditions in allConditions){
    cat(paste0("Analyzing condition(s): ",paste0(conditions,collapse=" vs. "),"\n"))
    #conditionFileIds = citrus.combinedFCSSet$fileIds[,conditions]
    
    # Calculate fold features
    #citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,minimumClusterSizePercent=minimumClusterSizePercent,medianColumns=medianColumns,mc.cores=4)
    #citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,minimumClusterSizePercent=minimumClusterSizePercent,conditions=conditions,medianColumns=medianColumns,mc.cores=4)
    cat("\tBuilding Fold Features\n")
    citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,conditions=conditions,...)
    conditionFeatures[[paste(rev(conditions),collapse="_vs_")]] = citrus.foldFeatureSet
    
    # Endpoint regress for each model type
    #citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
    #citrus.regressionResults = citrus.endpointRegress("glmnet",citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
    cat("\tAnalyzing vs. endpoint\n")
    citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family,...)
    names(citrus.regressionResults) = modelTypes
    conditionRegressionResults[[paste(rev(conditions),collapse="_vs_")]] = citrus.regressionResults
    
    # Plot results if requested
    if (plot){
      cat("\tPlotting Results\n")
      conditionOutputDir = file.path(outputDirectory,paste(rev(conditions),collapse="_vs_"))
      dir.create(conditionOutputDir,showWarnings=F)
      mclapply(citrus.regressionResults,citrus.plotRegressionResults,outputDirectory=conditionOutputDir,citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.combinedFCSSet=citrus.combinedFCSSet,family=family,labels=labels,conditions=conditions)
    }  
    cat("\n")
  }
  
  
  
  # Return Results
  results = list(citrus.foldClustering=citrus.foldClustering,conditionFoldFeatures=conditionFeatures,conditionRegressionResults=conditionRegressionResults)
  class(results) = "citrus.full.result"
  return(results)
}


#citrus.mapAndPredict = function(citrusResult,dataDir,newFileList,fileSampleSize,mappingColumns=NULL,transformCols=NULL,transformCofactor=5){
#  mappingResults = citrus.mapFileDataToClustering(dataDir=dataDir,newFileList=newFileList,fileSampleSize=fileSampleSize,preClusterResult=citrusResult$preClusterResult,mappingColumns=mappingColumns,transformCols=transformColumns,transformCofactor=transformCofactor)
#  modelLargeEnoughClusters = lapply(names(citrusResult),.extractConditionLargeEnoughClusters,foldFeatures=trainFeatures)
#  mappedFeatures = citrus.buildFeatures(preclusterResult=mappingResults,outputDir=outputDir,featureTypes=c("abundances","medians"),largeEnoughClusters=trainLargeEnoughClusters,medianColumns=medianCols)
#}