citrus.buildModel = function(features,labels,type,regularizationThresholds){
  if (type=="pamr"){
    model = pamr.train(list(x=t(features),y=labels),threshold=regularizationThresholds,remove.zeros=F)
  } else if (type=="glmnet") {
    # NOTE THAT THIS IS BINOMIAL EXPLICITLY. DOES MULTINOMIAL WORK THE SAME, IF ONLY 2 CLASSES PROVIDED?
    model = glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds)
  } else {
    stop(paste("Type:",type,"not yet implemented"));
  }
  return(model)
}

citrus.buildFoldModels = function(index,folds,foldFeatures,labels,type,regularizationThresholds){
  if (!((length(folds[[index]])==1) && (folds[[index]]=="all"))){
    labels = labels[-folds[[index]]]
  }
  citrus.buildModel(features=foldFeatures[[index]],labels=labels,type=type,regularizationThresholds=regularizationThresholds)
}

citrus.foldPredict = function(index,models,features){
  citrus.predict(models[[index]],features[[index]])
}

citrus.foldScore = function(index,folds,predictions,labels){
  return(predictions[[index]]==labels[folds[[index]]])
}

citrus.predict = function(model,features){
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


citrus.generateRgularizationThresholds = function(features,lables,modelTypes,n=50,alpha=1){
  if (length(modelTypes)<1){
    stop("no regularzation threshold types specified.")
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
    regs$glmnet = rev(glmnet(x=features,y=labels,family=family,alpha=alpha,nlambda=n)$lambda)
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


