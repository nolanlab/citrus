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
         regularizationThreshold=regularizationThresholds)
         
  class(foldModels) = "citrus.foldModels"
  return(foldModels)
}

citrus.endpointRegress = function(modelType,citrus.foldFeatureSet,labels,family,...){
  
  result = list()
  # Reg Thresholds
  result$regularizationThresholds = citrus.generateRegularizationThresholds.classification(features=citrus.foldFeatureSet$allFeatures,labels=labels,modelType=modelType,n=100)
  
  # Fold Models
  if (citrus.foldFeatureSet$nFolds>1){
    #foldModels = citrus.buildCrossValidatedEndpointModel(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=regularizationThresholds,family=family)
    result$foldModels = citrus.buildCrossValidatedEndpointModel(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=result$regularizationThresholds,family=family,...)
  } 
  
  # Final Models
  #result$finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=result$regularizationThresholds)
  result$finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=result$regularizationThresholds,...)
  
  # Calculate CV error rates
  if (citrus.foldFeatureSet$nFolds>1){
    result$thresholdCVRates = do.call(paste0("citrus.thresholdCVs.",family),args=list(foldModels=result$foldModels,leftoutFeatures=citrus.foldFeatureSet$leftoutFeatures,foldFeatures=citrus.foldFeatureSet$foldFeatures,modelType=modelType,regularizationThresholds=result$regularizationThresholds,labels=labels,folds=citrus.foldFeatureSet$folds))
  } else {
    result$thresholdCVRates = do.call(paste0("citrus.thresholdCVs.",family,".quick"),args=list(modelType=modelType,features=citrus.foldFeatureSet$allFeatures,labels=labels,regularizationThresholds=result$regularizationThresholds)) 
  }
  
  
  # Find CV Minima
  result$cvMinima = citrus.getCVMinima(modelType,thresholdCVRates=result$thresholdCVRates)
  
  # Extract differential features
  result$differentialFeatures = citrus.extractModelFeatures(cvMinima=result$cvMinima,finalModel=result$finalModel,finalFeatures=citrus.foldFeatureSet$allFeatures,regularizationThresholds=result$regularizationThresholds,family=family)
  
  # Extra info
  result$modelType=modelType
  result$family=family
  
  class(result) = "citrus.regressionResult"
  return(result)
}


citrus.thresholdCVs.classification.quick = function(modelType,features,labels,regularizationThresholds,nCVFolds=10,...){
  
  errorRates = list()
  errorRates$threshold=regularizationThresholds
  
  if (modelType=="pamr"){
    pamrData = list(x=t(features),y=labels)
    pamrModel = pamr.train(data=pamrData,threshold=regularizationThresholds,remove.zeros=F)
    pamrCVModel = pamr.cv(fit=pamrModel,data=pamrData,nfold=nCVFolds)
    errorRates$cvm = pamrCVModel$error
    
    cvmSD = as.vector(apply(sapply(pamrCVModel$folds,function(foldIndices,y,yhat){
      apply(apply(yhat[foldIndices,],2,"==",y[foldIndices]),2,sum)/length(foldIndices)
    },y=labels,yhat=pamrCVModel$yhat),1,sd))
    
    errorRates$cvsd = cvmSD / sqrt(length(pamrCVModel$folds))
    errorRates$fdr =  citrus:::pamr.fdr.new(pamrModel,data=pamrData,nperms=1000)$results[,"Median FDR"]
  } else if (modelType=="glmnet"){
    glmnetFamily="binomial"
    if (length(unique(labels))>2){
      glmnetFamily="multinomial"
    }
    addtlArgs = list(...)
    alpha=1
    if ("alpha" %in% names(addtlArgs)){
      alpha=addtlArgs[["alpha"]]
    }
    standardize=T
    if ("standardize" %in% names(addtlArgs)){
      standardize=addtlArgs[["standardize"]]
    }
    glmnetModel = cv.glmnet(x=features,y=labels,family=glmnetFamily,lambda=regularizationThresholds,type.measure="class",alpha=alpha,standardize=standardize)
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

###############################################
# Get the minima of models
###############################################
citrus.getCVMinima = function(modelType,thresholdCVRates,fdrRate=0.01){
  cvPoints=list();
  if (modelType=="sam"){
    cvPoints[["fdr_0.10"]]=10
    cvPoints[["fdr_0.05"]]=5
    cvPoints[["fdr_0.01"]]=1
  } else {
    errorRates = thresholdCVRates$cvm
    SEMs = thresholdCVRates$cvsd
    FDRRates = thresholdCVRates$fdr
    cvPoints[["cv.min.index"]] = min(which(errorRates==min(errorRates,na.rm=T)))
    cvPoints[["cv.min"]] = thresholdCVRates$threshold[cvPoints[["cv.min.index"]]]
    cvPoints[["cv.1se.index"]] = min(which(errorRates<=(errorRates[cvPoints[["cv.min.index"]]]+SEMs[cvPoints[["cv.min.index"]]])))
    cvPoints[["cv.1se"]] = thresholdCVRates$threshold[cvPoints[["cv.1se.index"]]]
    if (!is.null(FDRRates)) {
      if (any(FDRRates<fdrRate)){
        if (length(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates,na.rm=T))))>0){
          cvPoints[["cv.fdr.constrained.index"]] = max(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates,na.rm=T))))
          cvPoints[["cv.fdr.constrained"]] = thresholdCVRates$threshold[cvPoints[["cv.fdr.constrained.index"]]]
        }
      }
      
    }
  }
  return(cvPoints)
}

###############################################
# Get model features
###############################################
citrus.extractModelFeatures = function(cvMinima,finalModel,finalFeatures,regularizationThresholds,family){
  res = list();
  modelType = finalModel$type
  finalModel = finalModel$model
  for (cvPoint in names(cvMinima)[!grepl("index",names(cvMinima))]){
    threshold = cvMinima[[cvPoint]]
    thresholdIndex = cvMinima[[paste(cvPoint,"index",sep=".")]]
    if (modelType=="pamr"){
      if (finalModel$nonzero[thresholdIndex]>0){
        f = pamr.listgenes(fit=finalModel,data=list(x=t(finalFeatures),geneids=colnames(finalFeatures)),threshold=threshold)  
        f = as.vector(f[,1])
        res[[cvPoint]][["features"]] = f
        res[[cvPoint]][["clusters"]] = sort(unique(as.numeric(do.call("rbind",strsplit(f,split=" "))[,2])))  
      } else {
        res[[cvPoint]][["features"]] = NULL
        res[[cvPoint]][["clusters"]] = NULL
      }
      
    } else if (modelType=="glmnet"){
      # THIS NEEDS TO BE FIXED IN ORDER TO SUPPORT MULTINOMIAL REGRESSION WITH GLMNET
      f = as.matrix(predict(finalModel,newx=finalFeatures,type="coefficient",s=threshold))
      f = rownames(f)[f!=0]
      if ("(Intercept)" %in% f){
        f = f[-(which(f=="(Intercept)"))]
      }
      if (length(f)>0){
        res[[cvPoint]][["features"]] = f
        res[[cvPoint]][["clusters"]] = sort(unique(as.numeric(do.call("rbind",strsplit(f,split=" "))[,2])))  
      } else {
        res[[cvPoint]][["features"]] = NULL;
        res[[cvPoint]][["clusters"]] = NULL;
      }
    } else if (modelType=="sam"){
      sigGenes = rbind(finalModel$siggenes.table$genes.up,finalModel$siggenes.table$genes.lo)
      sigGenes = sigGenes[as.numeric(sigGenes[,"q-value(%)"])<threshold,,drop=F]
      f = sigGenes[,"Gene ID"]
      if (length(f)>0){
        #sigGenes = sigGenes[order(abs(as.numeric(sigGenes[,"Fold Change"]))),,drop=F]
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
