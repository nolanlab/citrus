citrus.buildModel.classification = function(features,labels,type,regularizationThresholds,...){
  
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
    model = pamr.train(data=pamrData,threshold=regularizationThresholds,remove.zeros=F)
  } else if (type=="glmnet") {
    # NOTE THAT THIS IS BINOMIAL EXPLICITLY. DOES MULTINOMIAL WORK THE SAME, IF ONLY 2 CLASSES PROVIDED?
    model = glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds,alpha=alpha,standardize=standardize)
  } else {
    stop(paste("Type:",type,"not yet implemented"));
  }
  return(model)
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

citrus.cvIteration.classification = function(i,modelType,features,labels,regularizationThresholds,nFolds,pamrModel=NULL,alpha=NULL,standardize=NULL){
  if (modelType == "pamr"){
    return(pamr.cv(fit=pamrModel,data=features,nfold=nFolds)$error)
  } else if (modelType=="glmnet"){
    return(cv.glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds,type.measure="class",nfolds=nFolds,alpha=alpha,standardize=standardize)$cvm)
  } else {
    stop(paste("Model Type",modelType,"unknown."));
  }
}


citrus.thresholdCVs.quick = function(foldModels,foldFeatures,modelTypes,regularizationThresholds,labels,family,...){
  sapply(modelTypes,citrus.thresholdCVs.model.quick,features=foldFeatures[[1]],regularizationThresholds=regularizationThresholds,family=family,labels=labels,...)
}

citrus.thresholdCVs.model.quick = function(modelType,features,regularizationThresholds,family,labels,nFolds=5,ncvIterations=10,...){
  typeRegularizationThresholds=regularizationThresholds[[modelType]]
  if (modelType=="pamr"){
    if (family=="survival"){
      stop("PAM model not implemeted for survival.")
    }
    pamrData = list(x=t(features),y=labels)
    pamrModel = pamr.train(data=pamrData,threshold=typeRegularizationThresholds,remove.zeros=F)
    errorRates = sapply(1:ncvIterations,paste("citrus.cvIteration",family,sep="."),modelType="pamr",features=pamrData,labels=labels,regularizationThresholds=typeRegularizationThresholds,nFolds=nFolds,pamrModel=pamrModel)
    fdrRate = citrus:::pamr.fdr.new(pamrModel,data=pamrData,nperms=1000)$results[,"Median FDR"]
  } else if (modelType=="glmnet"){
    addtlArgs = list(...)
    alpha=1
    if ("alpha" %in% names(addtlArgs)){
      alpha=addtlArgs[["alpha"]]
    }
    standardize=T
    if ("standardize" %in% names(addtlArgs)){
      standardize=addtlArgs[["standardize"]]
    }
    errorRates = sapply(1:ncvIterations,paste("citrus.cvIteration",family,sep="."),modelType="glmnet",features=features,labels=labels,regularizationThresholds=typeRegularizationThresholds,nFolds=nFolds,alpha=alpha,standardize=standardize)
  } else {
    stop(paste("CV for Model type",modelType,"not implemented"))
  }
  cvm=apply(errorRates,1,mean)
  cvsd=apply(errorRates,1,sd)/sqrt(ncvIterations)
  results = data.frame(threshold=typeRegularizationThresholds,cvm=cvm,cvsd=cvsd)
  if (exists("fdrRate")){
    results$fdr=fdrRate
  }
  return(results)
}

citrus.thresholdCVs.classification = function(foldModels,leftoutFeatures,foldFeatures,modelTypes,regularizationThresholds,labels,folds,...){
  leftoutPredictions = lapply(modelTypes,citrus.foldTypePredict,foldModels=foldModels,leftoutFeatures=leftoutFeatures)
  names(leftoutPredictions)=modelTypes
  
  predictionSuccess = lapply(as.list(modelTypes),citrus.foldTypeScore,folds=folds,leftoutPredictions=leftoutPredictions,labels=labels)
  names(predictionSuccess)=modelTypes
  
  thresholdErrorRates = lapply(modelTypes,citrus.calculateTypeErroRate,predictionSuccess=predictionSuccess,regularizationThresholds=regularizationThresholds)
  names(thresholdErrorRates)=modelTypes
  
  thresholdFDRRates = lapply(modelTypes,citrus.calculateTypeFDRRate,foldModels=foldModels,foldFeatures=foldFeatures,labels=labels)
  names(thresholdFDRRates)=modelTypes
  
  res=list()
  for (modelType in modelTypes){
    df = data.frame(threshold=regularizationThresholds[[modelType]],cvm=thresholdErrorRates[[modelType]]$cvm,cvsd=thresholdErrorRates[[modelType]]$cvsd);
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


citrus.generateRegularizationThresholds.classification = function(features,labels,modelTypes,n,...){
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
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
  regs = list()
  
  if ("pamr" %in% modelTypes){
    regs$pamr = rev(pamr.train(data=list(x=t(features),y=labels),n.threshold=n)$threshold)
  }
  if ("glmnet" %in% modelTypes){
    if (length(unique(labels))==2){
      glmfamily="binomial"
    } else {
      glmfamily="multinomial"
    }
    regs$glmnet = glmnet(x=features,y=labels,family=glmfamily,alpha=alpha,nlambda=n,standardize=standardize)$lambda
  }
  return(regs)
}

#citrus.calculateSEM = function(predictionSuccess,regularizationThresholds){
#  apply(do.call("rbind",lapply(predictionSuccess,citrus.getFoldErrorRate)),2,mean)/sqrt(nfolds)
#}

citrus.getFoldErrorRate = function(foldErrorRate){
  apply(!foldErrorRate,2,sum)/nrow(foldErrorRate)
}

citrus.buildModels=function(folds,foldFeatures,labels,regularizationThresholds,modelTypes,family,...){
  lapply(modelTypes,citrus.buildTypeModels,folds=folds,foldFeatures=foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,family=family,...=...)
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

#citrus.modelTypeSEM = function(modelType,predictionSuccess,regularizationThresholds){
#  citrus.calculateSEM(predictionSuccess[[modelType]],regularizationThresholds[[modelType]])
#}

citrus.calculateTypeErroRate = function(modelType,predictionSuccess,regularizationThresholds){
  nFolds=length(predictionSuccess[[modelType]])
  counter=1;
  tmp=list()
  for (i in 1:nFolds){
    for (j in 1:nrow(predictionSuccess[[modelType]][[i]])){
      tmp[[counter]] = predictionSuccess[[modelType]][[i]][j,]
      length(tmp[[counter]])=length(regularizationThresholds[[modelType]])
      counter=counter+1;
    }
  }
  bound = do.call("rbind",tmp)
  thresholdMeans= 1-apply(bound,2,mean,na.rm=T)
  thresholdSEMs = apply(bound,2,sd,na.rm=T)/sqrt(apply(!is.na(bound),2,sum))
  return(list(cvm=thresholdMeans,cvsd=thresholdSEMs))
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