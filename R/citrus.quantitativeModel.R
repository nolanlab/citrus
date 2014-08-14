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
  mse = sum((predictions[[index]]-labels[folds[[index]]])^2)/length(predictions[[index]])
  return(mse)
}

#' @rdname citrus.predict
#' @name citrus.predict
#' @export
citrus.predict.classification = function(citrus.endpointModel,newFeatures){
  if (citrus.endpointModel$type=="glmnet"){
    predictions = predict(citrus.endpointModel$model,newx=newFeatures,type="class")
  } else if (citrus.endpointModel$type=="pamr"){
    predictions = pamr.predictmany(fit=citrus.endpointModel$model,newx=t(newFeatures))$predclass
  } else {
    stop(paste("don't know how to predict for class",citrus.endpointModel$type));
  }
  rownames(predictions) = rownames(newFeatures)
  return(predictions)
}

#' @rdname citrus.generateRegularizationThresholds
#' @name citrus.generateRegularizationThresholds
#' @export
citrus.generateRegularizationThresholds.classification = function(features,labels,modelType,n=100,...){
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  
  if (modelType=="pamr"){
    return(rev(pamr.train(data=list(x=t(features),y=labels),n.threshold=n)$threshold))
  }
  
  if (modelType=="glmnet"){
    if (length(unique(labels))==2){
      glmfamily="binomial"
    } else {
      glmfamily="multinomial"
    }
    return(glmnet(x=features,y=labels,family=glmfamily,alpha=alpha,nlambda=n,standardize=standardize)$lambda)
  }
  
}


.calculatePredictionErrorRate = function(predictionSuccess,regularizationThresholds){
  nFolds=length(predictionSuccess)
  counter=1;
  tmp=list()
  for (i in 1:nFolds){
    for (j in 1:nrow(predictionSuccess[[i]])){
      tmp[[counter]] = predictionSuccess[[i]][j,]
      length(tmp[[counter]])=length(regularizationThresholds)
      counter=counter+1;
    }
  }
  bound = do.call("rbind",tmp)
  thresholdMeans= 1-apply(bound,2,mean,na.rm=T)
  thresholdSEMs = apply(bound,2,sd,na.rm=T)/sqrt(apply(!is.na(bound),2,sum))
  return(list(cvm=thresholdMeans,cvsd=thresholdSEMs))
} 

.calculateTypeFDRRate = function(foldModels,foldFeatures,labels,modelType){
  if (modelType=="pamr"){
    # FIX THIS
    # Should be average FDR across all models, not just all model 
    # return(pamr.fdr.new(foldModels[[modelType]][[length(foldModels[[modelType]])]],data=list(x=t(foldFeatures[[length(foldModels)]]),y=labels),nperms=1000)$results[,"Median FDR"])
    return(NULL)
  } else {
    return(NULL)
  }  
}
