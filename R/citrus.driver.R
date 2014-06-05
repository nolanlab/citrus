citrus.full = function(dataDirectory,
                       clusteringColumns,
                       labels,
                       family,
                       fileList,
                       outputDirectory,
                       modelTypes=c("glmnet"),
                       nFolds=1,
                       plot=T,
                       ...){
  
  # No point in running cv if SAM only model
  if (all(modelTypes=="sam")){
    nFolds=1
  }
  
  # Read in data
  citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory=dataDirectory,fileList=fileList,useChannelDescriptions=T,...)
    
  # Cluster each fold
  citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels=labels,nFolds=nFolds)
  save(list=c("citrus.combinedFCSSet","citrus.foldClustering"),file=file.path(outputDirectory,"citrusClustering.rData"),compress=F)
    
  # Calculate fold features
  #citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,minimumClusterSizePercent=minimumClusterSizePercent,medianColumns=medianColumns,mc.cores=4)
  citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,...)
  
  # Endpoint regress for each model type
  citrus.regressionResults = lapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
  names(citrus.regressionResults) = modelTypes
    
  # Plot results if requested
  if (plot){
    lapply(citrus.regressionResults,citrus.plotRegressionResults,outputDirectory=outputDirectory,citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.combinedFCSSet=citrus.combinedFCSSet,family=family,labels=labels)
  }
  
  # Return Results
  results = list(citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.regressionResults=citrus.regressionResults)
  class(results) = "citrus.full.result"
  return(results)
}


#citrus.mapAndPredict = function(citrusResult,dataDir,newFileList,fileSampleSize,mappingColumns=NULL,transformCols=NULL,transformCofactor=5){
#  mappingResults = citrus.mapFileDataToClustering(dataDir=dataDir,newFileList=newFileList,fileSampleSize=fileSampleSize,preClusterResult=citrusResult$preClusterResult,mappingColumns=mappingColumns,transformCols=transformColumns,transformCofactor=transformCofactor)
#  modelLargeEnoughClusters = lapply(names(citrusResult),.extractConditionLargeEnoughClusters,foldFeatures=trainFeatures)
#  mappedFeatures = citrus.buildFeatures(preclusterResult=mappingResults,outputDir=outputDir,featureTypes=c("abundances","medians"),largeEnoughClusters=trainLargeEnoughClusters,medianColumns=medianCols)
#}