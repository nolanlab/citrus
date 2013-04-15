citrus.generateRegularizationThresholds.survival = function(features,labels,modelTypes,n=100,alpha=1){
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
  }
  regs = list()
  if ("pamr" %in% modelTypes){
    stop("pamr model not implemented for survival data.");
  }
  if ("glmnet" %in% modelTypes){
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    regs$glmnet = rev(glmnet(x=features,y=s,family="cox",alpha=alpha,nlambda=c(n-1))$lambda)
    regs$glmnet[length(regs$glmnet)]=((regs$glmnet[length(regs$glmnet)-1]-regs$glmnet[length(regs$glmnet)-2])*1.5)+regs$glmnet[length(regs$glmnet)-1]
  }
  return(regs)
}

citrus.buildModel.survival = function(features,labels,type,regularizationThresholds,cv=F,nFolds=NULL,ncvRuns=10){
  
  if ((cv)&&(is.null(nFolds))){
    stop("nfolds not specififed for cross validation.")
  }
  
  if (type=="glmnet") {
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    glmmodel = glmnet(x=features,y=s,family="cox",lambda=regularizationThresholds,maxit=200000)
    if (cv){
      errorRates = sapply(1:ncvRuns,citrus.cvIteration.survival,modelType="glmnet",features=features,labels=labels,regularizationThresholds=regularizationThresholds,nFolds=nFolds)  
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

citrus.cvIteration.survival = function(i,modelType,features,labels,regularizationThresholds,nFolds){
  if (modelType=="glmnet"){
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    return(cv.glmnet(x=features,y=s,family="cox",lambda=regularizationThresholds,type.measure="class",nfolds=nFolds)$cvm)
  } else {
    stop(paste("Model Type",modelType,"unknown."));
  }
}

citrus.thresholdCVs.survival = function(foldModels,foldFeatures,leftoutFeatures,modelTypes,regularizationThresholds,labels,folds,...){
  
  s = Surv(time=labels[,"time"],event=labels[,"event"])
  
  thresholdPartialLikelihoods = lapply(modelTypes,citrus.calculateModelPartialLikelihood,leftoutFeatures=leftoutFeatures,foldFeatures=foldFeatures,foldModels=foldModels,labels=s,folds=folds,regularizationThresholds=regularizationThresholds)
  names(thresholdPartialLikelihoods)=modelTypes
  
  res=list()
  for (modelType in modelTypes){
    df = data.frame(threshold=regularizationThresholds[[modelType]],thresholdPartialLikelihoods[[modelType]]);
    res[[modelType]]=df
  }
  return(res)
}


citrus.calculateModelPartialLikelihood=function(modelType,leftoutFeatures,foldFeatures,foldModels,labels,folds,regularizationThresholds){
  nFolds=length(leftoutFeatures)
  nAllFolds = nFolds+1
  
  if (modelType=="glmnet"){
    cvRaw=matrix(NA,ncol=length(regularizationThresholds[[modelType]]),nrow=nFolds)
    for (foldId in 1:nFolds){
      coefmat = predict(foldModels[[modelType]][[foldId]], type = "coeff")
      plFull = coxnet.deviance(x=rbind(foldFeatures[[foldId]],leftoutFeatures[[foldId]]),y=rbind(labels[-folds[[foldId]],],labels[folds[[foldId]],]),offset=NULL,weights=rep(1,nrow(foldFeatures[[nAllFolds]])),beta=coefmat)
      plLeaveIn = coxnet.deviance(x=foldFeatures[[foldId]],y=labels[-folds[[foldId]],],offset=NULL,weights=rep(1,nrow(foldFeatures[[foldId]])),beta=coefmat)
      cvRaw[foldId,seq(along = plFull)]=plFull-plLeaveIn
    }
    
    foldid = rep(1,nrow(foldFeatures[[nAllFolds]]))
    for (i in 1:nFolds){
      foldid[folds[[i]]]=i
    }
    
    N = nFolds - apply(is.na(cvRaw), 2, sum)
    weights = as.vector(tapply(rep(1,nrow(foldFeatures[[nAllFolds]])) * labels[,"status"], foldid, sum))
    cvRaw = cvRaw/weights
    cvm=apply(cvRaw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvRaw, cvm, FALSE)^2, 2, weighted.mean, w = weights, na.rm = TRUE)/(N - 1))
    return(list(cvm=cvm,cvsd=cvsd))
  } else {
    stop(paste("Partial Likelihood calculations not implemented for model type",modelType));
  }
}

citrus.predict.survival = function(model,features){
  if ("glmnet" %in% class(model)){
    predictions = predict(model,newx=features)
  } else {
    stop(paste("don't know how to predict for class",class(model)));
  }
  return(predictions)
}