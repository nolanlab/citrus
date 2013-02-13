citrus.buildModel = function(features,labels,type,regularizationThresholds,cv=F,nFolds=NULL,ncvRuns=10){
  
  if ((cv)&&(is.null(nFolds))){
    stop("nfolds not specififed for cross validation.")
  }
  
  if (type=="pamr"){
    pamrData = list(x=t(features),y=labels)
    pamrModel = pamr.train(data=pamrData,threshold=regularizationThresholds,remove.zeros=F)
    if (cv){
      errorRates = sapply(1:ncvRuns,citrus.cvIteration,modelType="pamr",features=pamrData,labels=NULL,regularizationThresholds=regularizationThresholds,nFolds=nFolds,pamrModel=pamrModel)  
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
    glmmodel = glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds)
    if (cv){
      errorRates = sapply(1:ncvRuns,citrus.cvIteration,modelType="glmnet",features=features,labels=labels,regularizationThresholds=regularizationThresholds,nFolds=nFolds)  
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

citrus.cvIteration = function(i,modelType,features,labels,regularizationThresholds,nFolds,pamrModel=NULL){
  if (modelType == "pamr"){
    return(pamr.cv(fit=pamrModel,data=features,nfold=nFolds)$error)
  } else if (modelType=="glmnet"){
    return(cv.glmnet(x=features,y=labels,family="binomial",lambda=regularizationThresholds,type.measure="class",nfolds=nFolds)$cvm)
  } else {
    stop(paste("Model Type",modelType,"unknown."));
  }
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


citrus.generateRgularizationThresholds = function(features,labels,modelTypes,n=100,alpha=1){
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

citrus.buildTypeModels = function(modelType,folds,foldFeatures,labels,regularizationThresholds){
  lapply(1:length(folds),citrus.buildFoldModels,folds=folds,foldFeatures=foldFeatures,labels=labels,type=modelType,regularizationThresholds=regularizationThresholds[[modelType]])  
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

citrus.getCVMinima = function(modelType,thresholdErrorRates,thresholdSEMs,thresholdFDRRates){
  errorRates = thresholdErrorRates[[modelType]]
  SEMs = thresholdSEMs[[modelType]]
  FDRRates = thresholdFDRRates[[modelType]]
  cvPoints=list();
  cvPoints[["cv.min"]] = min(which(errorRates==min(errorRates)))
  cvPoints[["cv.1se"]] = min(which(errorRates<=(errorRates[cvPoints[["cv.min"]]]+SEMs[cvPoints[["cv.min"]]])))
  if (!is.null(FDRRates)) {
    if (any(FDRRates<0.01)){
      if (length(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates))))>0){
        cvPoints[["cv.fdr.constrained"]] = max(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates))))    
      }
    }
    
  }
  return(cvPoints)
}

citrus.extractModelFeatures = function(modelType,cvMinima,foldModels,foldFeatures,regularizationThresholds){
  res = list();
  nAllFolds = length(foldModels[[modelType]])
  finalModel = foldModels[[modelType]][[nAllFolds]]
  for (cvPoint in names(cvMinima[[modelType]])){
    if (modelType=="pamr"){
      threshold = regularizationThresholds[[modelType]][ cvMinima[[modelType]][[cvPoint]] ]
      if (finalModel$nonzero[cvMinima[[modelType]][[cvPoint]]]>0){
        f = pamr.listgenes(fit=finalModel,data=list(x=t(foldFeatures[[nAllFolds]]),geneids=colnames(foldFeatures[[nAllFolds]])),threshold=regularizationThresholds[[modelType]][ cvMinima[[modelType]][[cvPoint]] ] )  
        f = as.vector(f[,1])
        res[[cvPoint]][["features"]] = f
        res[[cvPoint]][["clusters"]] = sort(unique(as.numeric(do.call("rbind",strsplit(f,split=" "))[,2])))  
      } else {
        res[[cvPoint]][["features"]] = NULL
        res[[cvPoint]][["clusters"]] = NULL
      }
      
    } else if (modelType=="glmnet"){
      threshold = rev(regularizationThresholds[[modelType]])[ cvMinima[[modelType]][[cvPoint]] ]
      f = as.matrix(predict(finalModel,newx=foldFeatures[[nAllFolds]],type="coefficient",s=threshold))
      f = rownames(f)[f!=0][-1]
      if (length(f)>0){
        res[[cvPoint]][["features"]] = f
        res[[cvPoint]][["clusters"]] = sort(unique(as.numeric(do.call("rbind",strsplit(f,split=" "))[,2])))  
      } else {
        res[[cvPoint]][["features"]] = NULL;
        res[[cvPoint]][["clusters"]] = NULL;
      }
      
    }
  }
  return(res)
}

citrus.getModelTypes = function(){
  return(c("pamr","glmnet"))
}