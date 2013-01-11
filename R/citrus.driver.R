#' Driver to run a full citrus analysis
#' 
#' \code{citrus.full} runs a full Citrus analysis, including clustering, feature extraction and regularized classification. Plots are generated that show classification cross validation error rates, informative classifier features, and phenotype plots for informative clusters.
#' @param dataDir The path to the data directory containing the FCS files to be analyzed.
#' @param outputDir Path to a directory where the citrus output will be placed.
#' @param clusterCols A vector of integers or parameter names that should be used for clustering of the data
#' @param fileSampleSize Number of events to be sampled from each analyzed file. Files with fewer events contribute all of their events. 
#' @param fileList A matrix containing file names, conditions, and group assignments. See details. 
#' @param nFolds Number of cross validation folds used to assess model accuracy
#' @param modelTypes A vector of classification models to construct. Valid arguments are \code{pamr} and \code{glmnet}.
#' @param featureTypes A vector of descriptive feature types to be calculated for each cluster. Valid arguments are \code{densities} and \code{medians}. See details.
#' @param minimumClusterSizePercent Specifies the minimum cluster size to be analyzed as a percentage of the total aggregate datasize. A value etween \code{0} and \code{1}.
#' @param transformCols A vector of integer or parameter names to be transformed before analysis. 
#' @param ... Further arguments to be passed to Citrus subcomponents. 
#' @details Details about the cluster conditions matrix, fold features, etc.
#' @author Robert Bruggner
#' @references http://github.com/nolanlab/citrus/
citrus.full = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,transformCols=NULL,...){
  
  # Error check before we actually start the work.
  if ((!all(featureTypes %in% citrus.getFeatureTypes()))||(length(featureTypes)<1)){
    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.getFeatureTypes(),collapse=", "),"\n"))
  }
  
  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  
  labelCol = which(colnames(fileList)=="labels")
  if (length(labelCol)==0){
    stop("'labels' column not found in file list.")
  }
  
  if ("conditionComparaMatrix" %in% names(list(...))){
    #allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
    allConditions = citrus.convertConditionMatrix(list(...)[["conditionComparaMatrix"]]) 
  } else {
    allConditions = as.list(colnames(fileList)[-labelCol])
  }

  for (conditions in allConditions){
    cat(paste("Analyzing Condition",paste(conditions,collapse=" vs "),"\n"))
    
    conditionOutputDir = file.path(outputDir,paste(conditions,collapse="_vs_"))
    
    citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize)
    
    
    nAllFolds = nFolds+1
    folds = pamr:::balanced.folds(y=fileList[,labelCol],nfolds=nAllFolds)
    folds[[nAllFolds]]="all"
      
    
    cat("Clustering\n")
    foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=clusterCols,conditions=conditions)
    cat("Assigning Events to Clusters\n")
    foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)
    cat("Assigning Leftout Events to Clusters\n")
    leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterCols=clusterCols)
    cat("Calculating Fold Large Enough Clusters\n")
    foldLargeEnoughClusters = lapply(1:nAllFolds,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=foldsClusterAssignments,folds=folds,citrus.dataArray=citrus.dataArray,minimumClusterSizePercent=minimumClusterSizePercent)
      
    foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,...)
    #foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=medianColumns)
    leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,...)
    #leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,medianColumns=medianColumns)
    
    if (length(conditions)==2){
      foldFeatures = lapply(1:length(folds),citrus.buildFoldFeatureDifferences,features=foldFeatures,conditions=conditions,citrus.dataArray=citrus.dataArray,folds=folds)
      leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatureDifferences,features=leftoutFeatures,conditions=conditions,citrus.dataArray=citrus.dataArray,folds=folds,calculateLeaveoutData=T)
    }
    
    regularizationThresholds = citrus.generateRgularizationThresholds(foldFeatures[[nAllFolds]],labels=fileList[,labelCol],modelTypes=modelTypes)
    
    foldModels = lapply(modelTypes,citrus.buildTypeModels,folds=folds,foldFeatures=foldFeatures,labels=fileList[,labelCol],regularizationThresholds=regularizationThresholds)
    names(foldModels)=modelTypes
    
    leftoutPredictions = lapply(modelTypes,citrus.foldTypePredict,foldModels=foldModels,leftoutFeatures=leftoutFeatures)
    names(leftoutPredictions)=modelTypes
    
    predictionSuccess = lapply(as.list(modelTypes),citrus.foldTypeScore,folds=folds,leftoutPredictions=leftoutPredictions,labels=fileList[,labelCol])
    names(predictionSuccess)=modelTypes
    
    thresholdSEMs = lapply(modelTypes,citrus.modelTypeSEM,predictionSuccess=predictionSuccess)
    names(thresholdSEMs)=modelTypes
    
    thresholdErrorRates = lapply(modelTypes,citrus.calcualteTypeErroRate,predictionSuccess=predictionSuccess)
    names(thresholdErrorRates)=modelTypes
    
    thresholdFDRRates = lapply(modelTypes,citrus.calculateTypeFDRRate,foldModels=foldModels,foldFeatures=foldFeatures,labels=fileList[,labelCol])
    names(thresholdFDRRates)=modelTypes
    
    cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdErrorRates=thresholdErrorRates,thresholdSEMs=thresholdSEMs,thresholdFDRRates=thresholdFDRRates)
    names(cvMinima)=modelTypes
    
    # Make condition output directoy
    dir.create(conditionOutputDir,showWarnings=T,recursive=T)
    
    # Plot
    sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=regularizationThresholds,thresholdErrorRates=thresholdErrorRates,thresholdFDRRates=thresholdFDRRates,cvMinima=cvMinima,thresholdSEMs=thresholdSEMs,foldModels=foldModels)
    
    # Extract Features
    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
    #citrus.extractModelFeatures("pamr",cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
    names(differentialFeatures) = modelTypes
    
    # Plot Features
    lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=differentialFeatures,foldFeatures=foldFeatures,outputDir=conditionOutputDir,labels=fileList[,labelCol])
    
    # Plot Clusters
    lapply(modelTypes,citrus.plotClusters,differentialFeatures=differentialFeatures,outputDir=conditionOutputDir,clusterChildren=foldsClusterAssignments,citrus.dataArray=citrus.dataArray,conditions=conditions,clusterCols=clusterCols)
  }
}


citrus.quick = function(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds,modelTypes=c("pamr","glmnet"),featureTypes=c("densities"),minimumClusterSizePercent=0.05,transformCols=NULL,...){
  
  # Error check before we actually start the work.
  if ((!all(featureTypes %in% citrus.getFeatureTypes()))||(length(featureTypes)<1)){
    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.getFeatureTypes(),collapse=", "),"\n"))
  }
  
  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  if (!file.exists(outputDir)){
    stop(paste("Output directory",outputDir,"not found."))
  }
  
  labelCol = which(colnames(fileList)=="labels")
  if (length(labelCol)==0){
    stop("'labels' column not found in file list.")
  }
  
  if ("conditionComparaMatrix" %in% names(list(...))){
    #allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
    allConditions = citrus.convertConditionMatrix(list(...)[["conditionComparaMatrix"]]) 
  } else {
    allConditions = as.list(colnames(fileList)[-labelCol])
  }
  
  for (conditions in allConditions){
    cat(paste("Analyzing Condition",paste(conditions,collapse=" vs "),"\n"))
    
    conditionOutputDir = file.path(outputDir,paste(conditions,collapse="_vs_"))
    
    citrus.dataArray = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions=conditions,transformCols=transformCols,fileSampleSize=fileSampleSize)
    
    conditionData = citrus.dataArray$data[citrus.dataArray$data[,"fileId"]%in%citrus.dataArray$fileIds[,conditions],]
    cat("Clustering\n")
    clustering = citrus.cluster(data=conditionData[,clusterCols])
    cat("Assigning Events to Clusters\n")
    clusterAssignments = citrus.calculateCompleteHierarchicalMembership(clustering)
    cat("Calculating Large Enough Clusters\n")
    minimumClusterSize=sum(citrus.dataArray$data[,"fileId"]%in%citrus.dataArray$fileIds[,conditions])*minimumClusterSizePercent
    largeEnoughClusters = citrus.calculateLargeEnoughClusters(clusterAssignments,minimumClusterSize)

    features = citrus.buildFeatures(clusterAssignments=clusterAssignments,featureTypes=featureTypes,data=conditionData,largeEnoughClusters=largeEnoughClusters,foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions]),foldFileNames=citrus.dataArray$fileNames[as.vector(citrus.dataArray$fileIds[,conditions])],...)
    #features = citrus.buildFeatures(clusterAssignments=clusterAssignments,featureTypes=featureTypes,data=conditionData,largeEnoughClusters=largeEnoughClusters,foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions]),foldFileNames=citrus.dataArray$fileNames[as.vector(citrus.dataArray$fileIds[,conditions])],medianColumns=medianColumns)
    
    if (length(conditions)==2){
      features = features[citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[2]]],]-features[citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[1]]],]
      colnames(features) = paste(colnames(features),"difference")
    }
    
    regularizationThresholds = citrus.generateRgularizationThresholds(features,labels=fileList[,labelCol],modelTypes=modelTypes)
    
    models = list();
    thresholdErrorRates = list();
    thresholdSEMs = list();
    cvMinima = list();
    for (modelType in modelTypes){
      typeModel = citrus.buildModel(features=features,labels=fileList[,labelCol],type=modelType,regularizationThresholds=regularizationThresholds[[modelType]],cv=T,nFolds=nFolds)
      models[[modelType]] = list(typeModel$model)
      thresholdErrorRates[[modelType]]=typeModel$errorRates
      thresholdSEMs[[modelType]]= typeModel$se
      cvMinima[[modelType]]=list();
      cvMinima[[modelType]][["cv.min"]]=typeModel$cvmin
    }
        
    # Make condition output directoy
    dir.create(conditionOutputDir,showWarnings=T,recursive=T)
    
    # Plot
    sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=regularizationThresholds,thresholdErrorRates=thresholdErrorRates,thresholdFDRRates=NULL,cvMinima=cvMinima,thresholdSEMs=thresholdSEMs,foldModels=models)
    
    # Extract Features
    differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=models,foldFeatures=list(features),regularizationThresholds=regularizationThresholds)
    #citrus.extractModelFeatures("pamr",cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
    names(differentialFeatures) = modelTypes
    
    # Plot Features
    lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=differentialFeatures,foldFeatures=list(features),outputDir=conditionOutputDir,labels=fileList[,labelCol])
    
    # Plot Clusters
    lapply(modelTypes,citrus.plotClusters,differentialFeatures=differentialFeatures,outputDir=conditionOutputDir,clusterChildren=list(clusterAssignments),citrus.dataArray=citrus.dataArray,conditions=conditions,clusterCols=clusterCols)
  }
}


citrus.buildFoldFeatureDifferences = function(index,features,conditions,citrus.dataArray,folds,calculateLeaveoutData=F){
  
  if (folds[[index]]=="all"){
    includeRowIds=1:nrow(citrus.dataArray$fileIds)
  } else if (!calculateLeaveoutData){
    includeRowIds=setdiff(1:nrow(citrus.dataArray$fileIds),folds[[index]])
  } else {
    includeRowIds = folds[[index]]
  }
  diffFeatures = features[[index]][citrus.dataArray$fileNames[citrus.dataArray$fileIds[includeRowIds,conditions[2]]],] - features[[index]][citrus.dataArray$fileNames[citrus.dataArray$fileIds[includeRowIds,conditions[1]]],] 
  colnames(diffFeatures) = paste(colnames(diffFeatures),"difference")
  return(diffFeatures)
}
