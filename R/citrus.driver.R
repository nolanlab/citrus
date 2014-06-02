citrus.full = function(dataDirectory,
                       clusteringColumns,
                       labels,
                       family,
                       fileList,
                       outputDirectory,
                       modelTypes=c("glmnet"),
                       nFolds=1,
                       featureTypes=c("abundances"),
                       minimumClusterSizePercent=0.05,
                       fileSampleSize=NULL,
                       transformColumns=NULL,
                       conditionComparaMatrix=NULL,
                       plot=T,
                       transformCofactor=NULL,
                       ...){
  
  #balanceFactor=NULL
  #if (family=="survival"){
  #  balanceFactor=as.factor(labels[,"event"])
  #  if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
  #    stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
  #  }
  #}
  
  if (is.null(fileSampleSize)){
    fileSampleSize = 100/minimumClusterSizePercent
  }
  
  # No point in running cv if SAM only model
  if (all(modelTypes=="sam")){
    nFolds=1
  }
  
  # Read in data
  citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory=dataDirectory,fileList=fileList,fileSampleSize=fileSampleSize,transformColumns=transformColumns,useChannelDescriptions=T)
    
  # Cluster each fold
  citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels=labels,nFolds=nFolds)
  save(list=c("citrus.combinedFCSSet","citrus.foldClustering"),file=file.path(outputDirectory,"citrusClustering.rData"),compress=F)
    
  # Calculate fold features
  citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,minimumClusterSizePercent=minimumClusterSizePercent)
  
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