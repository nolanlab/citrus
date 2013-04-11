citrus.generateRegularizationThresholds.survival = function(features,labels,modelTypes,n=100,alpha=1){
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
  }
  regs = list()
  if ("pamr" %in% modelTypes){
    stop("pamr model not implemented for survival data.");
  }
  if ("glmnet" %in% modelTypes){
    s = Surv(time=labels[,1],event=labels[,2])
    regs$glmnet = rev(glmnet(x=features,y=s,family="cox",alpha=alpha,nlambda=n)$lambda)
  }
  return(regs)
}

citrus.buildModel.survival = function(features,labels,type,regularizationThresholds,cv=F,nFolds=NULL,ncvRuns=10){
  
  if ((cv)&&(is.null(nFolds))){
    stop("nfolds not specififed for cross validation.")
  }
  
  if (type=="glmnet") {
    s = Surv(time=labels[,"time"],event=labels[,"event"])
    glmmodel = glmnet(x=features,y=s,family="cox",lambda=regularizationThresholds)
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
