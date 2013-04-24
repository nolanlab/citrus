citrus.buildModel.classification = function(features,labels,type,regularizationThresholds,cv=F,nFolds=NULL,ncvRuns=10,...){
  
  if ((cv)&&(is.null(nFolds))){
    stop("nfolds not specififed for cross validation.")
  }
  
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  
  if (type=="pamr"){
    pamrData = list(x=t(features),y=labels)
    pamrModel = pamr.train(data=pamrData,threshold=regularizationThresholds,remove.zeros=F)
    if (cv){
      errorRates = sapply(1:ncvRuns,citrus.cvIteration.classification,modelType="pamr",features=pamrData,labels=NULL,regularizationThresholds=regularizationThresholds,nFolds=nFolds,pamrModel=pamrModel)  
      model=list();
      model$model=pamrModel;
      model$errorRates=apply(errorRates,1,mean)
      model$se = apply(errorRates,1,sd)/sqrt(ncvRuns)
      model$cvmin = min(which(model$errorRates==min(model$errorRates)))
    } else {
      model = pamrModel
    }
  } else if (type=="glmnet") {
    # NOTE THAT THIS IS BINOMIAL EXPLICITLY. DOES MULTINOMIAL WORK THE SAME, IF ONLY 2 CLASSES PROVIDED?
    glmmodel = glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds,alpha=alpha,standardize=standardize)
    if (cv){
      errorRates = sapply(1:ncvRuns,citrus.cvIteration.classification,modelType="glmnet",features=features,labels=labels,regularizationThresholds=regularizationThresholds,nFolds=nFolds,alpha=alpha,standardize=standardize)  
      model=list();
      model$model=glmmodel;
      model$errorRates=apply(errorRates,1,mean)
      model$se = apply(errorRates,1,sd)/sqrt(ncvRuns)
      model$cvmin = min(which(model$errorRates==min(model$errorRates)))
    } else {
      model = glmmodel
    }
  } else {
    stop(paste("Type:",type,"not yet implemented"));
  }
  return(model)
}

citrus.cvIteration.classification = function(i,modelType,features,labels,regularizationThresholds,nFolds,pamrModel=NULL,alpha=NULL,standardize=NULL){
  if (modelType == "pamr"){
    return(pamr.cv(fit=pamrModel,data=features,nfold=nFolds)$error)
  } else if (modelType=="glmnet"){
    return(cv.glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds,type.measure="class",nfolds=nFolds,alpha=alpha,standardize=standardize)$cvm)
  } else {
    stop(paste("Model Type",modelType,"unknown."));
  }
}

citrus.buildFoldModels = function(index,folds,foldFeatures,labels,type,regularizationThresholds,family,...){
  if (!((length(folds[[index]])==1) && (folds[[index]]=="all"))){
    if (!is.null(dim(labels))){
      labels = labels[-folds[[index]],]
    } else {
      labels = labels[-folds[[index]]]
    }
    
  }
  do.call(paste("citrus.buildModel",family,sep="."),args=list(features=foldFeatures[[index]],labels=labels,type=type,regularizationThresholds=regularizationThresholds,...=...))
}

citrus.thresholdCVs.classification = function(foldModels,leftoutFeatures,foldFeatures,modelTypes,regularizationThresholds,labels,folds,...){
  leftoutPredictions = lapply(modelTypes,citrus.foldTypePredict,foldModels=foldModels,leftoutFeatures=leftoutFeatures)
  names(leftoutPredictions)=modelTypes
  
  predictionSuccess = lapply(as.list(modelTypes),citrus.foldTypeScore,folds=folds,leftoutPredictions=leftoutPredictions,labels=labels)
  names(predictionSuccess)=modelTypes
  
  thresholdSEMs = lapply(modelTypes,citrus.modelTypeSEM,predictionSuccess=predictionSuccess)
  names(thresholdSEMs)=modelTypes
  
  thresholdErrorRates = lapply(modelTypes,citrus.calcualteTypeErroRate,predictionSuccess=predictionSuccess)
  names(thresholdErrorRates)=modelTypes
  
  thresholdFDRRates = lapply(modelTypes,citrus.calculateTypeFDRRate,foldModels=foldModels,foldFeatures=foldFeatures,labels=labels)
  names(thresholdFDRRates)=modelTypes
  
  res=list()
  for (modelType in modelTypes){
    df = data.frame(threshold=regularizationThresholds[[modelType]],cvm=thresholdErrorRates[[modelType]],cvsd=thresholdSEMs[[modelType]]);
    if (!is.null(thresholdFDRRates[[modelType]])){
      df = cbind(df,fdr=thresholdFDRRates[[modelType]]);  
    }
    res[[modelType]]=df
  }
  return(res)
}

citrus.foldPredict = function(index,models,features){
  citrus.predict.classification(models[[index]],features[[index]])
}

citrus.foldScore = function(index,folds,predictions,labels){
  return(predictions[[index]]==labels[folds[[index]]])
}

citrus.predict.classification = function(model,features){
  if ("glmnet" %in% class(model)){
    predictions = predict(model,newx=features,type="class")
  } else if (class(model)=="pamrtrained"){
    predictions = pamr.predictmany(fit=model,newx=t(features))$predclass
  } else {
    stop(paste("don't know how to predict for class",class(model)));
  }
  rownames(predictions) = rownames(features)
  return(predictions)
}


citrus.generateRegularizationThresholds.classification = function(features,labels,modelTypes,n=100,...){
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
  }
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  addtlArgs = list(...)
  if ((cv)&&(is.null(nFolds))){
    stop("nfolds not specififed for cross validation.")
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  regs = list()
  if ("pamr" %in% modelTypes){
    regs$pamr = rev(pamr.train(data=list(x=t(features),y=labels),n.threshold=n)$threshold)
  }
  if ("glmnet" %in% modelTypes){
    if (length(unique(labels))==2){
      family="binomial"
    } else {
      family="multinomial"
    }
    regs$glmnet = rev(glmnet(x=features,y=labels,family=family,alpha=alpha,nlambda=n,standardize=standardize)$lambda)
  }
  return(regs)
}

citrus.calculateSEM = function(predictionSuccess){
  nfolds = length(predictionSuccess)
  apply(do.call("rbind",lapply(predictionSuccess,citrus.getFoldErrorRate)),2,mean)/sqrt(nfolds)
}

citrus.getFoldErrorRate = function(foldErrorRate){
  apply(!foldErrorRate,2,sum)/nrow(foldErrorRate)
}

citrus.buildModels=function(folds,foldFeatures,labels,regularizationThresholds,modelTypes,family,...){
  lapply(modelTypes,citrus.buildTypeModels,folds=folds,foldFeatures=foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,family=family,...)
}

citrus.buildTypeModels = function(modelType,folds,foldFeatures,labels,regularizationThresholds,family,...){
  lapply(1:length(folds),citrus.buildFoldModels,folds=folds,foldFeatures=foldFeatures,labels=labels,type=modelType,regularizationThresholds=regularizationThresholds[[modelType]],family=family,...=...)
}

citrus.foldTypePredict = function(modelType,foldModels,leftoutFeatures){
  lapply(1:length(leftoutFeatures),citrus.foldPredict,models=foldModels[[modelType]],features=leftoutFeatures)
}
  
citrus.foldTypeScore = function(modelType,folds,leftoutPredictions,labels){
  lapply(1:length(leftoutPredictions[[modelType]]),citrus.foldScore,folds=folds,predictions=leftoutPredictions[[modelType]],labels=labels)
}

citrus.modelTypeSEM = function(modelType,predictionSuccess){
  citrus.calculateSEM(predictionSuccess[[modelType]])
}

citrus.calcualteTypeErroRate = function(modelType,predictionSuccess){
  return((1-(apply(do.call("rbind",predictionSuccess[[modelType]]),2,sum)/nrow(do.call("rbind",predictionSuccess[[modelType]])))  ))
} 

citrus.calculateTypeFDRRate = function(modelType,foldModels,foldFeatures,labels){
  if (modelType=="pamr"){
    return(pamr.fdr.new(foldModels[[modelType]][[length(foldModels[[modelType]])]],data=list(x=t(foldFeatures[[length(foldModels)]]),y=labels),nperms=1000)$results[,"Median FDR"])
  } else {
    return(NULL)
  }  
}

citrus.getModelTypes = function(){
  return(c("pamr","glmnet"))
}