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

citrus.buildFoldEndpointModel = function(foldIndex,folds,foldFeatures,labels,family,type,regularizationThreshold,...){
  foldLabels = labels[-folds[[foldIndex]]]
  citrus.buildEndpointModel(foldFeatures[[foldIndex]],labels=foldLabels,family=family,type=type,regularizationThreshold=regularizationThreshold,...)
}

citrus.buildCrossValidatedEndpointModel = function(type,citrus.foldFeatureSet,labels,regularizationThresholds,family="classification",...){
  addtlArgs = list(...)
    
  # Build models
  foldModels = lapply(1:citrus.foldFeatureSet$nFolds,
         citrus.buildFoldEndpointModel,
         folds=citrus.foldFeatureSet$folds,
         foldFeatures=citrus.foldFeatureSet$foldFeatures,
         labels=labels,
         family=family,
         type=type,
         regularizationThreshold=regularizationThresholds)
         
  class(foldModels) = "citrus.foldModels"
  return(foldModels)
}

citrus.endpointRegress = function(modelType,citrus.foldFeatureSet,labels,family,...){
  # Reg Thresholds
  regularizationThresholds = citrus.generateRegularizationThresholds.classification(features=citrus.foldFeatureSet$allFeatures,labels=labels,modelType=modelType,n=100)
  
  # Fold Models
  #foldModels = citrus.buildCrossValidatedEndpointModel(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=regularizationThresholds,family=family)
  foldModels = citrus.buildCrossValidatedEndpointModel(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=regularizationThresholds,family=family,...)
  
  # Final Models
  #finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=regularizationThresholds)
  finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=regularizationThresholds,...)
  
  # Calculate CV error rates
  thresholdCVRates = citrus.thresholdCVs.classification(foldModels=foldModels,leftoutFeatures=citrus.foldFeatureSet$leftoutFeatures,foldFeatures=citrus.foldFeatureSet$foldFeatures,modelType=modelType,regularizationThresholds=regularizationThresholds,labels=labels,folds=citrus.foldClustering$folds)
  
  # Find CV Minima
  cvMinima = citrus.getCVMinima(modelType,thresholdCVRates)
  
  # Extract differential features
  differentialFeatures = citrus.extractModelFeatures(cvMinima=cvMinima,finalModel=finalModel,finalFeatures=citrus.foldFeatureSet$allFeatures,regularizationThresholds=regularizationThresholds,family=family)
  
  result = list(differentialFeatures=differentialFeatures,cvMinima=cvMinima,thresholdCVRates=thresholdCVRates,foldModels=foldModels,finalModel=finalModel,regularizationThresholds=regularizationThresholds,modelType=modelType,family=family)
  class(result) = "citrus.regressionResult"
  return(result)
}

