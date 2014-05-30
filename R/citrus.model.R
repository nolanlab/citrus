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
         regularizationThreshold=regularizationThresholds[[type]])
         
  class(foldModels) = "citrus.foldModels"
  return(foldModels)
}

