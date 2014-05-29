citrus.buildEndpointModel = function(features,labels,family="classification",type="pamr",regularizationThresholds=NULL,...){
  if (is.null(regularizationThresholds)){
    regularizationThresholds = do.call(paste0("citrus.generateRegularizationThresholds.",family),args=list(features=features,labels=labels,modelType=type,n=100,...=...))
  }
  model = do.call(paste("citrus.buildModel",family,sep="."),args=list(features=features,labels=labels,type=type,regularizationThresholds=regularizationThresholds,...=...))
  result = list(model=model,regularizationThresholds=regularizationThresholds,family=family,type=type)
  class(result) = "citrus.endpointModel"
  return(result)  
}

print.citrus.endpointModel = function(x,...){
  cat("Citrus Model\n")
  cat(paste("\tFamily:",x$family,"\n"))
  cat(paste("\tType:",x$type,"\n"))
}

citrus.buildCrossValidatedEnpointModel(citrus.combinedFCSSet,clusteringColumns,labels,family,modelTypes,clusteringType="hierarchical",nFolds=10,...){
  addtlArgs = list(...)
  
  # Define Folds
  if ("nFolds" %in% names(addtlArgs)){
    nFolds = addtlArgs[["nFolds"]]
  }
  folds = pamr:::balanced.folds(y=labels,nfolds=nFolds)
  
  # Cluster fold events
  foldClustering = lapply(folds,citrus.clusterFold,citrus.combinedFCSSet=citrus.combinedFCSSet,clusteringColumns=clusteringColumns,clusteringType=clusteringType)
  allClustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns,clusteringType)
  
  # Map leftout data 
  foldMapping = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,folds=folds,foldClustering=foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,mc.cores=8)
  
  # Select clusters in each folds
  # Default is minimum cluster size
  foldLargeEnoughClusters = lapply(foldClustering,citrus.selectClusters,minimumClusterSizePercent=0.1)
  
  # Build Training Features
  foldFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,folds=folds,foldClusterIds=foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,foldClustering=foldClustering)
  
  # Build Testing Features
  
  # Build models
  
  # Predict testing data
  
  # Calculate error rates
  
  # Identify Critical CV Points
}

