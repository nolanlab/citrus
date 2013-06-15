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
  #featureObject = citrus.buildFeatures(preclusterResult,outputDir,featureTypes,minimumClusterSizePercent)
  featureObject = citrus.buildFeatures(preclusterResult,outputDir,featureTypes,minimumClusterSizePercent,...)
  #regressionResults = citrus.endpointRegress(citrus.featureObject=featureObject,family=family,modelTypes=modelTypes,labels=labels) 
  regressionResults = citrus.endpointRegress(citrus.featureObject=featureObject,family=family,modelTypes=modelTypes,labels=labels,...) 
  if (plot){
    citrus.plotRegressionResults(outputDir,citrus.preclusterResult=preclusterResult,citrus.featureObject=featureObject,citrus.regressionResult=regressionResults,modelTypes,family,labels)
  }
  return(list(preclusterResult=preclusterResult,featureObject=featureObject,regressionResults=regressionResults))
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
  
  regressionRes=list()
  for (conditionName in names(citrus.featureObject)){
    # Calculate Regularization Thresholds
    nAllFolds = length(citrus.featureObject[[conditionName]]$folds)
    if (lengthLabels!=nrow(citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]])){
      stop(paste("Number of Lables not equal to number of samples."))
    }    
    
    quick=F
    if (length(nAllFolds==1)&&citrus.featureObject[[conditionName]]$folds[[1]]=="all"){
      quick=T
    }
    
    cat("Calculating Regularization Thresholds\n")
    #regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]],labels=labels,modelTypes=modelTypes,n=100))
    regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]],labels=labels,modelTypes=modelTypes,n=100,...))
    
    # Build Fold Models
    cat("\nBuilding folds models\n")
    #foldModels = citrus.buildModels(folds=citrus.featureObject[[conditionName]]$folds,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family)
    foldModels = citrus.buildModels(folds=citrus.featureObject[[conditionName]]$folds,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family,...)
    names(foldModels)=modelTypes
    
    # Calculate fold deviances
    cat("\nCalculating threshold deviance rates\n")
    if (quick){
      thresholdCVRates = citrus.thresholdCVs.quick(foldModels=foldModels,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,modelTypes=modelTypes,regularizationThresholds=regularizationThresholds,labels=labels,family=family)
    } else {
      thresholdCVRates = do.call(paste("citrus.thresholdCVs",family,sep="."),args=list(foldModels=foldModels,leftoutFeatures=citrus.featureObject[[conditionName]]$leftoutFeatures,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,modelTypes=modelTypes,regularizationThresholds=regularizationThresholds,labels=labels,folds=citrus.featureObject[[conditionName]]$folds))  
    }
    
    # Find cv minima
    cat("Calculating CV minima\n")
    cvMinima = sapply(modelTypes,citrus.getCVMinima,thresholdCVRates=thresholdCVRates)
    
    # Extract Features
    cat("Extracting differential features\n")
    differentialFeatures = sapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=citrus.featureObject[[conditionName]]$foldFeatures,regularizationThresholds=regularizationThresholds,family=family)
    
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
