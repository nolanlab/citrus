#' Run a full citrus analysis
#' 
#' \code{citrus.full} runs a full Citrus analysis, including clustering, feature extraction and regularized classification. Plots are generated that show classification cross validation error rates, informative classifier features, and phenotype plots for informative clusters.
#' @param dataDir The path to the data directory containing the FCS files to be analyzed.
#' @param outputDir Path to a directory where the citrus output will be placed.
#' @param clusterCols A vector of integers or parameter names that should be used for clustering of the data
#' @param fileSampleSize Number of events to be sampled from each analyzed file. Files with fewer events contribute all of their events. 
#' @param labels Outcome variable to regress against. Should be groups for classification or survival time. 
#' @param nFolds Number of cross validation folds used to assess model accuracy. Set this value to "all" for quick analysis.
#' @param family Type of regression to perform. Valid options are \code{classification} and \code{survival}
#' @param fileList A matrix containing file names. 1 condition per column. See details.
#' @param filePopulationList A list with each named entry consisting of a matrix of file names. Samples in rows and 1 population per column.
#' @param modelTypes A vector of models to construct. Valid options are \code{pamr}, \code{glmnet},\code{sam}.
#' @param featureTypes A vector of descriptive feature types to be calculated for each cluster. Valid arguments are \code{abundances}, \code{medians}, and \code{emDists}. See details.
#' @param minimumClusterSizePercent Specifies the minimum cluster size to be analyzed as a percentage of the total aggregate datasize. A value etween \code{0} and \code{1}.
#' @param transformCols A vector of integer or parameter names to be transformed before analysis. 
#' @param plot Logical value indicating whether or not Citrus output plots should be created. Defaults to \code{TRUE}.
#' @param ... Further arguments to be passed to Citrus subcomponents. 
#' @details Details about the cluster conditions matrix, fold features, etc.
#' @author Robert Bruggner
#' @references http://github.com/nolanlab/citrus/
citrus.full = function(dataDirectory,
                       clusterCols,
                       labels,
                       family,
                       fileList,
                       outputDirectory,
                       modelTypes=c("glmnet"),
                       nFolds=1,
                       featureTypes=c("abundances"),
                       minimumClusterSizePercent=0.05,
                       fileSampleSize=NULL,
                       transformColumns=NULL,
                       conditionComparaMatrix=NULL,
                       plot=T,
                       transformCofactor=NULL,
                       ...){
  
  #balanceFactor=NULL
  #if (family=="survival"){
  #  balanceFactor=as.factor(labels[,"event"])
  #  if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
  #    stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
  #  }
  #}
  
  if (is.null(fileSampleSize)){
    fileSampleSize = 100/minimumClusterSizePercent
  }
  
  # No point in running cv if SAM only model
  if (all(modelTypes=="sam")){
    nFolds=1
  }
  
  # Read in data
  citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory=dataDirectory,fileList=fileList,fileSampleSize=fileSampleSize,transformColumns=transformColumns,useChannelDescriptions=T)
    
  # Cluster each fold
  citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels=labels,nFolds=nFolds)
  save(list=c("citrus.combinedFCSSet","citrus.foldClustering"),file=file.path(outputDirectory,"citrusClustering.rData"),compress=F)
    
  # Calculate fold features
  citrus.foldFeatureSet = citrus.buildFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,minimumClusterSizePercent=minimumClusterSizePercent)
  
  citrus.regressionResults = lapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
  names(citrus.regressionResults) = modelTypes
    
      
  if (plot){
    lapply(citrus.regressionResults,citrus.plotRegressionResults,outputDirectory=outputDirectory,citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.combinedFCSSet=citrus.combinedFCSSet,family=family,labels=labels)
  }
  results = list(citrus.foldClustering=citrus.foldClustering,citrus.foldFeatureSet=citrus.foldFeatureSet,citrus.regressionResults=citrus.regressionResults)
  class(results) = "citrus.full.result"
}

#citrus.endpointRegress = function(conditionFeatureList,family,modelTypes,labels,...){
#  
#  if (!(family %in% citrus.familyList())){
#    stop("'family' argument must specified and one of the following: ",paste(citrus.familyList(),collapse=", "))
#  }
#  
#  if (family=="survival" && ("pamr" %in% modelTypes)){
#    warning("'pamr' model not implemented for 'survival' analysis. Removing.");
#    modelTypes=modelTypes[-which(modelTypes=="pamr")]
#  }
#  
#  if (length(modelTypes)==0){
#    stop("No model types specified.")
#  }
#  
#  
#  if (family=="survival"){
#    if ((ncol(labels)!=2)||(!all(colnames(labels) %in% c("time","event")))){
#      stop("Incorrect labeling for files. Expecting 'time' and 'event' label columns.")
#    }
#    lengthLabels = nrow(labels)
#  } else {
#    lengthLabels = length(labels)
#  }
#  
#  regressionRes=list()
#  for (conditionName in names(conditionFeatureList)){
    # Calculate Regularization Thresholds
#    nAllFolds = length(conditionFeatureList[[conditionName]]$folds)
#    if (lengthLabels!=nrow(conditionFeatureList[[conditionName]]$foldFeatures[[nAllFolds]])){
#      stop(paste("Number of Lables not equal to number of samples."))
#    }    
    
#    quick=F
#    if (length(nAllFolds==1)&&conditionFeatureList[[conditionName]]$folds[[1]]=="all"){
#      quick=T
#    }
#    cat("Calculating Regularization Thresholds\n")
#    #regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=conditionFeatureList[[conditionName]]$foldFeatures[[nAllFolds]],labels=labels,modelTypes=modelTypes,n=100))
#    regularizationThresholds = do.call(paste("citrus.generateRegularizationThresholds",family,sep="."),args=list(features=conditionFeatureList[[conditionName]]$foldFeatures[[nAllFolds]],labels=labels,modelTypes=modelTypes,n=100,...))
#    
#    # Build Fold Models
#    cat("\nBuilding folds models\n")
#    #foldModels = citrus.buildModels(folds=conditionFeatureList[[conditionName]]$folds,foldFeatures=conditionFeatureList[[conditionName]]$foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family)
#    foldModels = citrus.buildModels(folds=conditionFeatureList[[conditionName]]$folds,foldFeatures=conditionFeatureList[[conditionName]]$foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,modelTypes=modelTypes,family=family,...)
#    names(foldModels)=modelTypes
#    
#    # Calculate fold deviances
#    cat("\nCalculating threshold deviance rates\n")
#    if (quick){
#      thresholdCVRates = citrus.thresholdCVs.quick(foldModels=foldModels,foldFeatures=conditionFeatureList[[conditionName]]$foldFeatures,modelTypes=modelTypes,regularizationThresholds=regularizationThresholds,labels=labels,family=family)
#    } else {
#      thresholdCVRates = do.call(paste("citrus.thresholdCVs",family,sep="."),args=list(foldModels=foldModels,leftoutFeatures=conditionFeatureList[[conditionName]]$leftoutFeatures,foldFeatures=conditionFeatureList[[conditionName]]$foldFeatures,modelTypes=modelTypes,regularizationThresholds=regularizationThresholds,labels=labels,folds=conditionFeatureList[[conditionName]]$folds))  
#    }
    
    # Find cv minima
#    cat("Calculating CV minima\n")
#    cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdCVRates=thresholdCVRates)
#    names(cvMinima) = modelTypes
#    
#    # Extract Features
#    cat("Extracting differential features\n")
#    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=conditionFeatureList[[conditionName]]$foldFeatures,regularizationThresholds=regularizationThresholds,family=family)
#    names(differentialFeatures) = modelTypes
#    
#    regressionRes[[conditionName]] = list(differentialFeatures=differentialFeatures,cvMinima=cvMinima,thresholdCVRates=thresholdCVRates,foldModels=foldModels,regularizationThresholds=regularizationThresholds)
#    
#  }
#  return(regressionRes) 
#}

#citrus.mapAndPredict = function(citrusResult,dataDir,newFileList,fileSampleSize,mappingColumns=NULL,transformCols=NULL,transformCofactor=5){
#  mappingResults = citrus.mapFileDataToClustering(dataDir=dataDir,newFileList=newFileList,fileSampleSize=fileSampleSize,preClusterResult=citrusResult$preClusterResult,mappingColumns=mappingColumns,transformCols=transformColumns,transformCofactor=transformCofactor)
#  modelLargeEnoughClusters = lapply(names(citrusResult),.extractConditionLargeEnoughClusters,foldFeatures=trainFeatures)
#  mappedFeatures = citrus.buildFeatures(preclusterResult=mappingResults,outputDir=outputDir,featureTypes=c("abundances","medians"),largeEnoughClusters=trainLargeEnoughClusters,medianColumns=medianCols)
#}