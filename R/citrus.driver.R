#' Run a full citrus analysis
#' 
#' \code{citrus.full} runs a full Citrus analysis, including clustering, feature extraction and regularized classification. Plots are generated that show classification cross validation error rates, informative classifier features, and phenotype plots for informative clusters.
#' @param dataDir The path to the data directory containing the FCS files to be analyzed.
#' @param outputDir Path to a directory where the citrus output will be placed.
#' @param clusterCols A vector of integers or parameter names that should be used for clustering of the data
#' @param fileSampleSize Number of events to be sampled from each analyzed file. Files with fewer events contribute all of their events. 
#' @param fileList A matrix containing file names, conditions, and response variables. See details
#' @param nFolds Number of cross validation folds used to assess model accuracy
#' @param family Specifies the type of response variable is being regressed. Valid options are \code{classification} or \code{survival}. See details for required parameters in file list.
#' @param modelTypes A vector of classification models to construct. Valid options are \code{pamr} and \code{glmnet}.
#' @param featureTypes A vector of descriptive feature types to be calculated for each cluster. Valid arguments are \code{densities} and \code{medians}. See details.
#' @param minimumClusterSizePercent Specifies the minimum cluster size to be analyzed as a percentage of the total aggregate datasize. A value etween \code{0} and \code{1}.
#' @param transformCols A vector of integer or parameter names to be transformed before analysis. 
#' @param plot Logical value indicating whether or not Citrus output plots should be created. Defaults to \code{TRUE}.
#' @param returnResults Logical value indicating whether or not caluclated clusters and features should be returned. Defaults to FALSE.
#' @param ... Further arguments to be passed to Citrus subcomponents. 
#' @details Details about the cluster conditions matrix, fold features, etc.
#' @author Robert Bruggner
#' @references http://github.com/nolanlab/citrus/
citrus.full = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,nFolds,family,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,transformCols=NULL,conditionComparaMatrix=NULL,plot=T,returnResults=F,transformFactor=NULL,...){
  balanceFactor=NULL
  if (family=="survival"){
    balanceFactor=as.factor(labels[,"event"])
    if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
      stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
    }
  }
  preclusterResult = citrus.preCluster(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,transformCols,conditionComparaMatrix,balanceFactor,transformFactor)
  featureObject = citrus.buildFeatures(preclusterResult,outputDir,featureTypes,minimumClusterSizePercent,...)
  regressionResults = citrus.endpointRegress(citrus.featureObject=featureObject,citrus.preclusterResult=preclusterResult,family,modelTypes,labels,...) 
  if (plot){
    citrus.plotRegressionResults(outputDir,citrus.preclusterResult=preclusterResult,citrus.featureObject=featureObject,citrus.regressionResult=regressionResults,modelTypes,family,labels)
  }
  return(list(preclusterResult=preclusterResult,featureObject=featureObject,regressionResults=regressionResults))
}

citrus.quick = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,family,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,transformCols=NULL,conditionComparaMatrix=NULL,plot=T,returnResults=F,transformFactor=NULL,...){
  balanceFactor=NULL
  if (family=="survival"){
    balanceFactor=as.factor(labels[,"event"])
    if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
      stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
    }
  }
  preclusterResult = citrus.preCluster(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds="all",transformCols,conditionComparaMatrix,balanceFactor,transformFactor)
  regresssionResult = citrus.endpointRegress.quick(preclusterResult,outputDir,family,labels,modelTypes,featureTypes,minimumClusterSizePercent,plot,returnResults,plot=plot,returnResults=returnResults,...)
  if (returnResults){
    return(regresssionResult)
  }
}


citrus.preCluster = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,transformCols=NULL,conditionComparaMatrix=NULL,balanceFactor=NULL,transformFactor=5){
  
  res=list()
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
  } else {
    allConditions = as.list(colnames(fileList))
  }
  
  if (is.null(balanceFactor)){
    balanceFactor = as.factor(sample(c(0,1),nrow(fileList),replace=T))
  } else if (levels(balanceFactor)==1){
    balanceFactor = as.factor(sample(c(0,1),length(balanceFactor),replace=T))
  }
  if (nFolds=="all"){
    folds = list()
    nAllFolds=1
  } else {
    folds = pamr:::balanced.folds(y=balanceFactor,nfolds=nFolds)
    nAllFolds = nFolds+1  
  }
  folds[[nAllFolds]]="all"
  
  for (conditions in allConditions){
    cat(paste("Clustering Condition",paste(conditions,collapse=" vs "),"\n"))
    
    citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize,transformFactor=transformFactor)
    
    foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions)
    cat("Assigning Events to Clusters\n")
    foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)
    resObj = list(folds=folds,foldsCluster=foldsCluster,foldsClusterAssignments=foldsClusterAssignments,conditions=conditions,citrus.dataArray=citrus.dataArray,clusterColumns=clusterCols)
    if (nFolds!="all"){
      cat("Assigning Leftout Events to Clusters\n")
      leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols)
      resObj$leftoutClusterAssignments=leftoutClusterAssignments
    }
    save(resObj,file=file.path(outputDir,paste("citrus.Cluster.",paste(conditions,collapse="_vs_"),".rDat",sep="")))
    res[[paste(conditions,collapse="_vs_")]]=resObj
  }
  return(res)
}

citrus.endpointRegress.quick = function(preclusterResult,outputDir,family,labels,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,plot=T,returnResults=F,...){
  
  addtlArgs = list(...)
  
  # Error check before we actually start the work.
  if ((!all(featureTypes %in% citrus.featureTypes()))||(length(featureTypes)<1)){
    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
  }
  
  if (!(family %in% citrus.familyList())){
    stop("'family' argument must specified and one of the following: ",paste(citrus.familyList(),collapse=", "))
  }
  
  if (family=="survival" && ("pamr" %in% modelTypes)){
    warning("'pamr' model not implemented for 'survival' analysis. Removing.");
    modelTypes=modelTypes[-which(modelTypes=="pamr")]
  }
  
  if (length(modelTypes)==0){
    stop("No model types specified.")
  }
  
  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  if (family=="survival"){
    if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
      stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
    }
    lengthLabels = nrow(labels)
  } else {
    lengthLabels = length(labels)
  }
  
  regressionRes = list()
  
  for (conditionName in names(preclusterResult)){
    
    if (lengthLabels!=nrow(preclusterResult[[conditionName]]$citrus.dataArray$fileIds)){
      stop(paste("Number of Lables not equal to number of samples:",lengthLabels,nrow(preclusterResult[[conditionName]]$citrus.dataArray$fileIds)))
    }    
    
    cat(paste("Analyzing Condition",conditionName,"\n"))
    conditions=preclusterResult[[conditionName]]$conditions
    
    cat("Calculating Large Enough Clusters\n")
    largeEnoughClusters = citrus.calculateFoldLargeEnoughClusters(1,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,folds=list("all"),citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,minimumClusterSizePercent=minimumClusterSizePercent)
    
    cat("Calculating Features\n")
    #features = citrus.buildFoldFeatures(1,featureTypes=featureTypes,folds=list("all"),citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=list(largeEnoughClusters),conditions=conditions,emdColumns=medianCols,medianColumns=medianCols)
    features = citrus.buildFoldFeatures(1,featureTypes=featureTypes,folds=list("all"),citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=list(largeEnoughClusters),conditions=conditions,...)
    if (is.null(dim(features))){
      stop("No Features Calculated.")
    }
    
    # Calculate Regularization Thresholds
    cat("Calculating Regularization Thresholds\n")
    #regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=features,labels=labels,modelTypes=modelTypes,n=100))
    regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=features,labels=labels,modelTypes=modelTypes,n=100,...))

    # Build Fold Models
    cat("\nBuilding models\n")
    #models = citrus.buildModels(folds=list("all"),foldFeatures=list(features),labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family)
    models = citrus.buildModels(folds=list("all"),foldFeatures=list(features),labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family,...)
    names(models)=modelTypes
      
    # Calculate fold deviances
    cat("\nCalculating threshold deviance rates\n")
    thresholdCVRates = lapply(modelTypes,citrus.thresholdCVs.quick,features=features,regularizationThresholds=regularizationThresholds,family=family,labels=labels)
    names(thresholdCVRates) = modelTypes
    
    # Find cv minima
    cat("Calculating CV minima\n")
    cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdCVRates=thresholdCVRates)
    names(cvMinima)=modelTypes
    
    # Extract Features
    cat("Extracting differential features\n")
    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=models,foldFeatures=list(features),regularizationThresholds=regularizationThresholds,family=family)
    names(differentialFeatures) = modelTypes
    differentialFeatures
    
    if (returnResults){
      regressionRes[[conditionName]] = list(largeEnoughClusters=largeEnoughClusters,features=features,differentialFeatures=differentialFeatures,cvMinima=cvMinima,thresholdCVRates=thresholdCVRates,models=models,regularizationThresholds=regularizationThresholds)
    }
    
    # Plot
    if (plot){
      # Make condition output directoy
      conditionOutputDir = file.path(outputDir,conditionName)
      dir.create(conditionOutputDir,showWarnings=T,recursive=T)
      
      # Plot
      cat("Plotting Model Error Rates\n")
      sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=regularizationThresholds,thresholdCVRates=thresholdCVRates,cvMinima=cvMinima,foldModels=models,family=family)
      
      # Plot Features
      cat("Plotting Stratifying Features\n")
      lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=differentialFeatures,features=features,outputDir=conditionOutputDir,labels=labels,family=family,cvMinima=cvMinima,foldModels=models,regularizationThresholds=regularizationThresholds)
      
      # Plot Clusters
      cat("Plotting Stratifying Clusters\n")
      lapply(modelTypes,citrus.plotClusters,differentialFeatures=differentialFeatures,outputDir=conditionOutputDir,clusterChildren=preclusterResult[[conditionName]]$foldsClusterAssignments,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,conditions=preclusterResult[[conditionName]]$conditions,clusterCols=preclusterResult[[conditionName]]$clusterColumns)
      
      # GO BACK AND MAKE THIS PARALLEL CALLS....
      cat("Plotting Clustering Graph\n")
      nAllFolds=1
      g = citrus.createHierarchyGraph(largeEnoughClusters=largeEnoughClusters,mergeOrder=preclusterResult[[conditionName]]$foldsCluster[[nAllFolds]]$merge,clusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]])
      l = layout.reingold.tilford(g,root=length(largeEnoughClusters),circular=T)
      clusterMedians = t(sapply(largeEnoughClusters,.getClusterMedians,clusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]],data=preclusterResult[[conditionName]]$citrus.dataArray$data,clusterCols=preclusterResult[[conditionName]]$clusterColumns))
      rownames(clusterMedians) = largeEnoughClusters
      for (modelType in names(differentialFeatures)){
        citrus.plotHierarchicalClusterMedians(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),"markerPlots.pdf"),clusterMedians,graph=g,layout=l)  
        for (cvPoint in names(differentialFeatures[[modelType]])){
          featureClusterMatrix = .getClusterFeatureMatrix(differentialFeatures[[modelType]][[cvPoint]][["features"]])
          citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),paste("featurePlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,largeEnoughClusters=largeEnoughClusters,graph=g,layout=l)
          citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),paste("featurePetalPlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,largeEnoughClusters=largeEnoughClusters,graph=g,layout=l,petalPlot=T,clusterMedians=clusterMedians)
        }
      }
    }
  }
  return(regressionRes) 
}


citrus.buildFeatures = function(preclusterResult,outputDir,featureTypes=c("densities"),minimumClusterSizePercent=0.05,returnResults=T,...){  
  addtlArgs = list(...)
  
  # Error check before we actually start the work.
  if ((!all(featureTypes %in% citrus.featureTypes()))||(length(featureTypes)<1)){
    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
  }
  
  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  
  
  featureRes = list()
  
  for (conditionName in names(preclusterResult)){
    
    cat(paste("Building features for condition",conditionName,"\n"))
    conditions=preclusterResult[[conditionName]]$conditions
    
    folds = preclusterResult[[conditionName]]$folds
    nAllFolds=length(folds)
    nFolds=nAllFolds-1
    
    cat("Calculating Fold Large Enough Clusters\n")
    foldLargeEnoughClusters = lapply(1:nAllFolds,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,minimumClusterSizePercent=minimumClusterSizePercent)
    
    cat("Calculating Features\n")
    foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,...)
    if (any(do.call("c",lapply(lapply(foldFeatures,dim),is.null)))){
      stop("No Features Calculated.")
    }
    #foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,emdColumns=emdColumns)
    leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,...)
    #leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,emdColumns=emdColumns)
    
    #Normalize features... Sometimes EMD's aren't calculated. Need a better way to handle this.
    nof = lapply(1:nFolds,citrus.getNonOverlappingFeatures,foldFeatures=foldFeatures,leftoutFeatures=leftoutFeatures)
    foldFeatures[1:nFolds] = lapply(1:nFolds,citrus.removeFeatures,foldFeatures=foldFeatures,nonOverlappingFeatures=nof)
    leftoutFeatures = lapply(1:nFolds,citrus.removeFeatures,foldFeatures=leftoutFeatures,nonOverlappingFeatures=nof)
    
    conditionResults = list(foldLargeEnoughClusters=foldLargeEnoughClusters,foldFeatures=foldFeatures,leftoutFeatures=leftoutFeatures,folds=folds)
    save(conditionResults,file=file.path(outputDir,paste("citrus.Features.",conditionName,".rDat",sep="")))
    featureRes[[conditionName]]=conditionResults
  }
  
  if (returnResults){
    return(featureRes)
  } else {
    return()
  }
  
}
  

citrus.endpointRegress = function(citrus.featureObject,family,modelTypes,labels,returnResults=T,...){
  
  if (!(family %in% citrus.familyList())){
    stop("'family' argument must specified and one of the following: ",paste(citrus.familyList(),collapse=", "))
  }
  
  if (family=="survival" && ("pamr" %in% modelTypes)){
    warning("'pamr' model not implemented for 'survival' analysis. Removing.");
    modelTypes=modelTypes[-which(modelTypes=="pamr")]
  }
  
  if (length(modelTypes)==0){
    stop("No model types specified.")
  }
  
  
  if (family=="survival"){
    if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
      stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
    }
    lengthLabels = nrow(labels)
  } else {
    lengthLabels = length(labels)
  }
  
  
  regressionRes = list()

  regressionRes=list()
  for (conditionName in names(citrus.featureObject)){
    # Calculate Regularization Thresholds
    nAllFolds = length(citrus.featureObject[[conditionName]]$foldFeatures)
    if (lengthLabels!=nrow(citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]])){
      stop(paste("Number of Lables not equal to number of samples."))
    }    
    
    cat("Calculating Regularization Thresholds\n")
    regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]],labels=labels,modelTypes=modelTypes,n=100,...))
    
    # Build Fold Models
    cat("\nBuilding folds models\n")
    foldModels = citrus.buildModels(folds=citrus.featureObject[[conditionName]]$folds,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family,...)
    names(foldModels)=modelTypes
    
    # Calculate fold deviances
    cat("\nCalculating threshold deviance rates\n")
    thresholdCVRates = do.call(paste("citrus.thresholdCVs",family,sep="."),args=list(foldModels=foldModels,leftoutFeatures=citrus.featureObject[[conditionName]]$leftoutFeatures,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,modelTypes=modelTypes,regularizationThresholds=regularizationThresholds,labels=labels,folds=citrus.featureObject[[conditionName]]$folds))
    
    # Find cv minima
    cat("Calculating CV minima\n")
    cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdCVRates=thresholdCVRates)
    names(cvMinima)=modelTypes
    
    # Extract Features
    cat("Extracting differential features\n")
    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,regularizationThresholds=regularizationThresholds,family=family)
    names(differentialFeatures) = modelTypes
  
    regressionRes[[conditionName]] = list(differentialFeatures=differentialFeatures,cvMinima=cvMinima,thresholdCVRates=thresholdCVRates,foldModels=foldModels,regularizationThresholds=regularizationThresholds)
    
  }
  if (returnResults){
    return(regressionRes)
  } else {
    return()
  }
}


# Plot
citrus.plotRegressionResults = function(outputDir,citrus.preclusterResult,citrus.featureObject,citrus.regressionResult,modelTypes,family,labels,...){
  addtlArgs = list(...)
  clusterColLabels=NULL;
  if ("clusterColLabels" %in% names(addtlArgs)){
    clusterColLabels=addtlArgs[["clusterColLabels"]]
  }
  for (conditionName in names(citrus.regressionResult)){
    nAllFolds = length(citrus.featureObject[[conditionName]]$foldFeatures)
    # Make condition output directoy
    conditionOutputDir = file.path(outputDir,conditionName)
    dir.create(conditionOutputDir,showWarnings=T,recursive=T)
  
    cat("Plotting Error Rate\n")
    sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds,thresholdCVRates=citrus.regressionResult[[conditionName]]$thresholdCVRates,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,family=family)
   
    cat("Plotting Stratifying Features\n")
    lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,features=citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]],outputDir=conditionOutputDir,labels=labels,family=family,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds)
    
    cat("Plotting Stratifying Clusters\n")
    lapply(modelTypes,citrus.plotClusters,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,outputDir=conditionOutputDir,clusterChildren=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments,citrus.dataArray=citrus.preclusterResult[[conditionName]]$citrus.dataArray,conditions=citrus.preclusterResult[[conditionName]]$conditions,clusterCols=citrus.preclusterResult[[conditionName]]$clusterColumns,clusterColLabels=clusterColLabels)
    
    
    # GO BACK AND MAKE THIS PARALLEL CALLS....
    cat("Plotting Clustering Graph\n")
    g = citrus.createHierarchyGraph(largeEnoughClusters=citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],mergeOrder=citrus.preclusterResult[[conditionName]]$foldsCluster[[nAllFolds]]$merge,clusterAssignments=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]])
    l = layout.reingold.tilford(g,root=length(citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]]),circular=T)
    clusterMedians = t(sapply(citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],.getClusterMedians,clusterAssignments=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]],data=citrus.preclusterResult[[conditionName]]$citrus.dataArray$data,clusterCols=citrus.preclusterResult[[conditionName]]$clusterColumns))
    rownames(clusterMedians) = citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]]
    for (modelType in names(citrus.regressionResult[[conditionName]]$differentialFeatures)){
      citrus.plotHierarchicalClusterMedians(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),"markerPlots.pdf"),clusterMedians,graph=g,layout=l)  
      for (cvPoint in names(citrus.regressionResult[[conditionName]]$differentialFeatures[[modelType]])){
        featureClusterMatrix = .getClusterFeatureMatrix(citrus.regressionResult[[conditionName]]$differentialFeatures[[modelType]][[cvPoint]][["features"]])
        citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),paste("featurePlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,largeEnoughClusters=citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],graph=g,layout=l)
        citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),paste("featurePetalPlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,largeEnoughClusters=citrus.featureObject[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],graph=g,layout=l,petalPlot=T,clusterMedians=clusterMedians)
      }
    }
  }
  
}
