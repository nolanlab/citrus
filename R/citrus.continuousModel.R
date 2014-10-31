#' @rdname citrus.buildEndpointModel
#' @name citrus.buildEndpointModel
#' @export
citrus.buildModel.continuous = function(features,labels,type,regularizationThresholds,...){
  
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  
  if (all(c("thisFoldIndex","finalModelIndex") %in% names(addtlArgs))){
    if ((type=="sam")&&(addtlArgs[["thisFoldIndex"]]!=addtlArgs[["finalModelIndex"]])){
      return(NULL)
    }
  }
  
  if (type=="glmnet") {
    model = glmnet(x=features,y=labels,family="gaussian",lambda=regularizationThresholds,alpha=alpha,standardize=standardize)
  } else if (type=="sam"){
    noVarianceFeatures = apply(features,2,var)==0
    model = SAM(x=t(features[,!noVarianceFeatures]),y=labels,resp.type="Quantitative",genenames=colnames(features[,!noVarianceFeatures]),nperms=10000)
  } else {
    stop(paste("Type:",type,"not implemented for continuous model"));
  }
  return(model)
}

#' @rdname citrus.thresholdCVs
#' @name citrus.thresholdCVs
#' @export
citrus.thresholdCVs.quick.continuous = function(modelType,features,labels,regularizationThresholds,nCVFolds=10,...){
  
  errorRates = list()
  errorRates$threshold=regularizationThresholds
  
  if (modelType=="glmnet"){
    addtlArgs = list(...)
    alpha=1
    if ("alpha" %in% names(addtlArgs)){
      alpha=addtlArgs[["alpha"]]
    }
    standardize=T
    if ("standardize" %in% names(addtlArgs)){
      standardize=addtlArgs[["standardize"]]
    }
    glmnetModel = cv.glmnet(x=features,y=labels,family="gaussian",lambda=regularizationThresholds,type.measure="mse",alpha=alpha,standardize=standardize)
    errorRates$cvm = glmnetModel$cvm
    errorRates$cvsd = glmnetModel$cvsd
  } else if (modelType=="sam"){
    warning("No thresholds for SAM. This is Normal.")
    return(NA)
  } else {
    stop(paste("CV for Model type",modelType,"not implemented"))
  }
  
  return(data.frame(errorRates))
}

foldPredict.continuous = function(index,models,features){
  citrus.predict.continuous(models[[index]],features[[index]])
}

foldScore.continuous = function(index,folds,predictions,labels){
  squaredDifference = (predictions[[index]]-labels[folds[[index]]])^2
  return(squaredDifference)
}

#' @rdname citrus.predict
#' @name citrus.predict
#' @export
citrus.predict.continuous = function(citrus.endpointModel,newFeatures){
  if (citrus.endpointModel$type=="glmnet"){
    predictions = predict(citrus.endpointModel$model,newx=newFeatures,type="response")
  } else {
    stop(paste("don't know how to predict for class",citrus.endpointModel$type));
  }
  rownames(predictions) = rownames(newFeatures)
  return(predictions)
}

#' @rdname citrus.generateRegularizationThresholds
#' @name citrus.generateRegularizationThresholds
#' @export
citrus.generateRegularizationThresholds.continuous = function(features,labels,modelType,n=100,...){
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  
  if (modelType=="glmnet"){
    return(glmnet(x=features,y=labels,family="gaussian",alpha=alpha,nlambda=n,standardize=standardize)$lambda)
  } else {
    return(NULL)
  }
}




