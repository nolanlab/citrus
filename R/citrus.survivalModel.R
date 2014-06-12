citrus.generateRegularizationThresholds.survival = function(features,labels,modelTypes,n,...){
  stop("Survival currently unsupported.")
  addtlArgs = list(...)
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
  }
  regs = list()
  if ("pamr" %in% modelTypes){
    stop("pamr model not implemented for survival data.");
  }
  if ("glmnet" %in% modelTypes){
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    regs$glmnet = glmnet(x=features,y=s,family="cox",alpha=alpha,nlambda=c(n-1),standardize=standardize)$lambda
    regs$glmnet[1]=((regs$glmnet[1]-regs$glmnet[2])*1.5)+regs$glmnet[1]
  }
  return(regs)
}

citrus.buildModel.survival = function(features,labels,type,regularizationThresholds,...){
  stop("Survival currently unsupported.")
  addtlArgs = list(...)
  alpha=1
  if ("alpha" %in% names(addtlArgs)){
    alpha = addtlArgs[["alpha"]]
  }
  standardize=T
  if ("standardize" %in% names(addtlArgs)){
    standardize=addtlArgs[["standardize"]]
  }
  if (type=="glmnet") {
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    model = glmnet(x=features,y=s,family="cox",lambda=regularizationThresholds,maxit=200000,alpha=alpha,standardize=standardize)
  } else {
    stop(paste("Type:",type,"not yet implemented"));
  }
  return(model)
}

citrus.cvIteration.survival = function(i,modelType,features,labels,regularizationThresholds,nFolds,alpha,standardize){
  if (modelType=="glmnet"){
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    return(cv.glmnet(x=features,y=s,family="cox",lambda=regularizationThresholds,nfolds=nFolds,alpha=alpha,standardize=standardize)$cvm)
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

citrus.predict.survival = function(model,features,s){
  if ("glmnet" %in% class(model)){
    predictions = predict(model,newx=features,s=s)
  } else {
    stop(paste("don't know how to predict for class",class(model)));
  }
  return(predictions)
}