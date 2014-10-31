#' Build an endpoint model 
#' 
#' This function constructs an endpoint model using features calculated by citrus. 
#' @name citrus.buildEndpointModel
#' @param features A numeric matrix of predictive features. Rows are observations and column entries are features. 
#' @param labels A vector of endpoint values (i.e. class labels) for each row of the feature matrix. 
#' @param family Family of endpoint model to be constructed. Valid values are \code{classification} and \code{continuous}. 
#' @param type Statistical model to be used. For \code{family="classification"}, options are \code{pamr} (Nearest Shrunken Centroid), \code{glmnet} (Lasso-regularized logistic regression), and \code{sam} (Non-parametric test in differences of means). For \code{family="continuous"}, options are \code{glmnet} (L1-regularized linear regression), and \code{sam}. 
#' @param regularizationThresholds Vector of regularization values for penalized model construction. If \code{NULL}, values are automatically generated. Not valid for \code{sam} models.
#' @param ... Other parameters passed to model-fitting procedures. 
#' 
#' @return An object of class \code{citrus.endpointModel} with properties:
#' \item{model}{The statistical model fit on supplied data.}
#' \item{regularizationThresholds}{Regularization Thresholds used to constrain penalized models.}
#' \item{family}{Family of model.}
#' \item{type}{Model type.}
#' 
#' @author Robert Bruggner
#' @export
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Build model
#' endpointModel = citrus.buildEndpointModel(abundanceFeatures,labels)
citrus.buildEndpointModel = function(features,labels,family="classification",type="pamr",regularizationThresholds=NULL,...){
  if (is.null(regularizationThresholds)){
    regularizationThresholds = citrus.generateRegularizationThresholds(features=features,labels=labels,modelType=type,family=family,...)
  }
  model = do.call(paste("citrus.buildModel",family,sep="."),args=list(features=features,labels=labels,type=type,regularizationThresholds=regularizationThresholds,...=...))
  result = list(model=model,regularizationThresholds=regularizationThresholds,family=family,type=type)
  class(result) = "citrus.endpointModel"
  return(result)  
}

#' @export
#' @name citrus.buildEndpointModel
print.citrus.endpointModel = function(citrus.endpointModel,...){
  cat("Citrus Model\n")
  cat(paste("\tFamily:",citrus.endpointModel$family,"\n"))
  cat(paste("\tType:",citrus.endpointModel$type,"\n"))
}

#' Generate model regularization thresholds 
#' 
#' Generate a range of regularization thresholds for model construction
#' @param features Features used to construct model
#' @param labels Endpoint lables for samples and features
#' @param modelType Method used to construct endpoint model. Valid options are: \code{pamr} and \code{glmnet}.
#' @param family Model family. Valid options are: \code{classification} and \code{continuous}.
#' @param n Number of regularization thresholds to generate
#' @param ... Other arguments passed to model-fitting methods
#' 
#' @return A vector of regularization threshold values.
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Calculate regularization thresholds
#' regularizationThresholds = citrus.generateRegularizationThresholds(abundanceFeatures,labels,modelType="pamr",family="classification")
citrus.generateRegularizationThresholds = function(features,labels,modelType,family,n=100,...){
  do.call(paste0("citrus.generateRegularizationThresholds.",family),args=list(features=features,labels=labels,modelType=modelType,n=n,...=...))
}

#' Calculate model error rates 
#' 
#' Calculate model error rates at different regularization thresholds. 
#' 
#' @name citrus.thresholdCVs
#' @param modelType Type of model to be constructed. Valid options are: \code{pamr} and \code{glmnet}.
#' @param foldFeatures List of features with each entry containing features from an independent clustering.
#' @param features Features calculated from a clustering of all samples.
#' @param labels Endpoint labels of clustered samples.
#' @param regularizationThresholds Thresholds for model regularization.
#' @param family Model family. Valid options are \code{classification} and \code{continuous}.
#' @param folds List of fold indices
#' @param foldModels Models constructed from each fold of features.
#' @param leftoutFeatures Features calculated for leftout samples mapped to clustered data space.
#' @param nCVFolds Number of folds for quick cross-validation.
#' @param ... Other parameters passsed to model-fitting methods.
#' 
#' @details If independent fold-clustering and fold-features are calculated, use \code{citrus.thresholdCVs}. 
#' If features are derived from a clustering of all samples together, use \code{citrus.thresholdCVs.quick}. See examples.
#' 
#' @return Matrix of model error rates, standard error of error estimates, and false discovery rates (if possible) at 
#' supplied regularization thresholds.
#' 
#' @author Robert Bruggner
#' @export 
#' 
#' @examples
#' ########################################
#' # Example of citrus.thresholdCVs.quick
#' ########################################
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Calculate regularization thresholds
#' regularizationThresholds = citrus.generateRegularizationThresholds.classification(abundanceFeatures,labels,modelType="pamr")
#' 
#' # Calculate CV Error rates
#' thresholdCVRates = citrus.thresholdCVs.quick("pamr",abundanceFeatures,labels,regularizationThresholds,family="classification") 
#' 
#' ########################################
#' # Example of citrus.thresholdCVs
#' ########################################
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster each fold
#' citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels,nFolds=4)
#' 
#' # Build fold features and leftout features
#' citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet)
#' 
#' # Build fold models 
#' citrus.foldModels = citrus.buildFoldsEndpointModels(type="pamr",citrus.foldFeatureSet,labels)
#' 
#' citrus.thresholdCVs(modelType="pamr",
#'                     foldFeatures=citrus.foldFeatureSet$foldFeatures,
#'                     labels=labels,
#'                     regularizationThresholds=citrus.foldModels[[1]]$regularizationThresholds,
#'                     family="classification",
#'                     folds=citrus.foldFeatureSet$folds,
#'                     foldModels=citrus.foldModels,
#'                     leftoutFeatures=citrus.foldFeatureSet$leftoutFeatures)
citrus.thresholdCVs = function(modelType,foldFeatures,labels,regularizationThresholds,family,folds,foldModels,leftoutFeatures,...){
  if (modelType=="sam"){
    return(NULL)
  }
  #do.call(paste0("citrus.thresholdCVs.",family),args=list(modelType=modelType,foldFeatures=foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds,folds=folds,foldModels=foldModels,leftoutFeatures=leftoutFeatures,...=...))
  leftoutPredictions = lapply(1:length(leftoutFeatures),paste0("foldPredict.",family),models=foldModels,features=leftoutFeatures)
  predictionScore = lapply(1:length(leftoutPredictions),paste0("foldScore.",family),folds=folds,predictions=leftoutPredictions,labels=labels)
  thresholdErrorRates = calculatePredictionErrorRate(predictionScore,regularizationThresholds,family)
  thresholdFDRRates = .calculateTypeFDRRate(foldModels=foldModels,foldFeatures=foldFeatures,labels=labels,modelType=modelType)  
  results = data.frame(threshold=regularizationThresholds,cvm=thresholdErrorRates$cvm,cvsd=thresholdErrorRates$cvsd);
  if (!is.null(thresholdFDRRates)){
    results$fdr = thresholdFDRRates
  }
  return(results)
}


#' @rdname citrus.thresholdCVs
#' @export 
citrus.thresholdCVs.quick = function(modelType,features,labels,regularizationThresholds,family,nCVFolds=10,...){
  if (modelType=="sam"){
    return(NULL)
  }
  do.call(paste0("citrus.thresholdCVs.quick.",family),args=list(modelType=modelType,features=features,labels=labels,regularizationThresholds=regularizationThresholds,nCVFolds=nCVFolds,...=...))  
}

#' Predict labels of new feature set
#' 
#' Predict labels of new feature set
#' @param citrus.endpointModel A \code{citrus.endpointModel} object.
#' @param newFeatures Features from samples to predict labels for. 
#' 
#' @return Matrix of predicted sample endpoints at all model regularization thresholds.
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # List of files to be clustered
#' fileList1 = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs")[seq(from=2,to=20,by=2)])
#' 
#' # List of files to be mapped
#' fileList2 = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs")[seq(from=1,to=19,by=2)])
#' 
#' # Read the data 
#' citrus.combinedFCSSet1 = citrus.readFCSSet(dataDirectory,fileList1)
#' citrus.combinedFCSSet2 = citrus.readFCSSet(dataDirectory,fileList2)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster first dataset
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet1,clusteringColumns)
#' 
#' # Map new data to exsting clustering
#' citrus.mapping = citrus.mapToClusterSpace(citrus.combinedFCSSet.new=citrus.combinedFCSSet2,citrus.combinedFCSSet.old=citrus.combinedFCSSet1,citrus.clustering)
#' 
#' # Large Enough Clusters 
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Clustered Features and mapped features
#' clusteredFeatures = citrus.calculateFeatures(citrus.combinedFCSSet1,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' mappedFeatures = citrus.calculateFeatures(citrus.combinedFCSSet2,clusterAssignments=citrus.mapping$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # Labels
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Build Endpoint Model 
#' citrus.endpointModel = citrus.buildEndpointModel(clusteredFeatures,labels[seq(from=2,to=20,by=2)])
#' 
#' # Predict
#' citrus.predict(citrus.endpointModel,newFeatures=mappedFeatures)
citrus.predict = function(citrus.endpointModel,newFeatures){
  do.call(paste0("citrus.predict.",citrus.endpointModel$family),args=list(citrus.endpointModel=citrus.endpointModel,newFeatures=newFeatures))
}
  
#' Build models from each fold of clustering
#' 
#' Builds a model from features derived from each independent fold of clustering.
#' 
#' @param type Model Type. Valid options are \code{pamr}, \code{glmnet}, and \code{sam}.
#' @param citrus.foldFeatureSet. A \code{citrus.foldFeatureSet} object.
#' @param labels Endpoint labels for samples. 
#' @param regularizationThresholds Regularization thresholds for penalized models. 
#' @param family Family of model to be constructed. Valid options are \code{classification} and \code{continuous}.
#' @param ... Other arguments passed to model-fitting functions.
#' 
#' @return A list of models, one model fit on each fold's feature set. 
#' 
#' @author Robert Bruggner
#' @export
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster each fold
#' citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels,nFolds=4)
#' 
#' # Build fold features and leftout features
#' citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet)
#' 
#' # Build fold models 
#' citrus.foldModels = citrus.buildFoldsEndpointModels(type="pamr",citrus.foldFeatureSet,labels)
citrus.buildFoldsEndpointModels = function(type,citrus.foldFeatureSet,labels,regularizationThresholds=NULL,family="classification",...){
  
  if (is.null(regularizationThresholds)){
    regularizationThresholds = citrus.generateRegularizationThresholds(features=citrus.foldFeatureSet$allFeatures,labels=labels,modelType=type,family=family,...)
  }
    
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

citrus.buildFoldEndpointModel = function(foldIndex,folds,foldFeatures,labels,family,type,regularizationThreshold,...){
  foldLabels = labels[-folds[[foldIndex]]]
  citrus.buildEndpointModel(foldFeatures[[foldIndex]],labels=foldLabels,family=family,type=type,regularizationThreshold=regularizationThreshold,...)
}


#' Regress against an experimental endpoint
#'
#' Regress cluster properties against an experimental endpoint of interest. Models are fit on supplied features and constrained 
#' by regularization thresholds (\code{glmnet} and \code{pamr}) or FDR (\code{sam}). Stratifying features are returned along with 
#' corresponding cluster IDs. 
#' 
#' @param modelType Method to be used for model-fitting. Valid options are: \code{glmnet},\code{pamr}, and \code{sam}.
#' @param citrus.foldFeatureSet A \code{citrus.foldFeatureSet} object.
#' @param labels Vector of endpoint values for analyzed samples.
#' @param family Family of model to fit. Valid options are \code{classification} and \code{continuous}.
#' @param ... Other parameters passed to model-fitting methods.
#' 
#' @details If independent clusterings are run (i.e. \code{citrus.clusterAndMapFolds} is run with \code{nFolds > 1}), model are fit on each 
#' feature set calculated for each clustering fold and final regularization thresholds are selected by predicting endpoint values for leftout samples whose data
#' was mapped to existing cluster space. If a single clustering was run (i.e. \code{citrus.clusterAndMapFolds} is run with \code{nFolds = 1}), 
#' cross-validation is used to select final regularization thresholds based on features derived from a clustering of all samples. Regardless
#' of how regularization thresholds are selected, the final reported features are from the final model constructed from all features, constrained by 
#' identified optimal regularization thresholds.
#' 
#' @return A \code{citrus.regression} object with the following properties:
#' \item{regularizationThresholds}{Regularization thresholds used to constrain all constructed models. Not applicable for \code{sam} models.}
#' \item{foldModels}{A \code{citrus.endpointModel} constructed from each independent fold feature set. \code{NULL} if \code{nFolds = 1}.}
#' \item{finalModel}{A \code{citrus.endpointModel} constructed from features derived from the clustering of all samples together.}
#' \item{thresholdCVRates}{Matrix containing the average error rate and standard error of models at each regularization threshold. FDR also reported where possible.}
#' \item{cvMinima}{Values and indicies of pre-selected cross-validation error-rate thresholds.}
#' \item{differentialFeatures}{Non-zero model features and corresponding clusters from the \code{finalModel} constrained by \code{cvMinima}.}
#' \item{modelType}{Type of model fit on data.}
#' \item{family}{Family of regression model.}
#' \item{labels}{Endpoint labels of analyzed samples.}
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Cluster data
#' citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels,nFolds=4)
#' 
#' # Build abundance features
#' citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet)
#' 
#' # Endpoint regress
#' citrus.regressionResult = citrus.endpointRegress(modelType="pamr",citrus.foldFeatureSet,labels,family="classification")
citrus.endpointRegress = function(modelType,citrus.foldFeatureSet,labels,family,...){
    
  if (nrow(citrus.foldFeatureSet$allFeatures)!=length(labels)){
    stop(paste0("Number of features (",nrow(citrus.foldFeatureSet$allFeatures),") different from length of labels (",length(labels),")."))
  }
    
  # Build results 
  result = list()
  
  # Reg Thresholds
  result$regularizationThresholds = citrus.generateRegularizationThresholds(features=citrus.foldFeatureSet$allFeatures,labels=labels,modelType=modelType,family=family,...)
  
  # Fold Models
  if ((citrus.foldFeatureSet$nFolds>1)&&(modelType!="sam")){
    #foldModels = citrus.buildFoldsEndpointModels(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=regularizationThresholds,family=family)
    result$foldModels = citrus.buildFoldsEndpointModels(type=modelType,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,regularizationThresholds=result$regularizationThresholds,family=family,...)
  } 
  
  # Final Models
  #result$finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=result$regularizationThresholds)
  result$finalModel = citrus.buildEndpointModel(features=citrus.foldFeatureSet$allFeatures,labels=labels,family=family,type=modelType,regularizationThresholds=result$regularizationThresholds,...)
  
  # Calculate CV error rates
  if (citrus.foldFeatureSet$nFolds>1){
    result$thresholdCVRates = citrus.thresholdCVs(modelType=modelType,
                                                  foldFeatures=citrus.foldFeatureSet$foldFeatures,
                                                  labels=labels,
                                                  regularizationThresholds=result$regularizationThresholds,
                                                  family=family,
                                                  folds=citrus.foldFeatureSet$folds,
                                                  foldModels=result$foldModels,
                                                  leftoutFeatures=citrus.foldFeatureSet$leftoutFeatures)
  } else {
    result$thresholdCVRates = citrus.thresholdCVs.quick(modelType=modelType,
                                                        features=citrus.foldFeatureSet$allFeatures,
                                                        labels=labels,
                                                        regularizationThresholds=result$regularizationThresholds,
                                                        family=family) 
  }
  
  
  # Find CV Minima
  result$cvMinima = citrus.getCVMinima(modelType,thresholdCVRates=result$thresholdCVRates)
  
  # Extract differential features
  result$differentialFeatures = citrus.extractModelFeatures(cvMinima=result$cvMinima,finalModel=result$finalModel,finalFeatures=citrus.foldFeatureSet$allFeatures)
  
  # Extra info
  result$modelType=modelType
  result$family=family
  result$labels=labels
  
  class(result) = "citrus.regressionResult"
  return(result)
}


#' Get regularization thresholds of pre-selected cross-validation points
#' 
#' #' Get regularization thresholds of pre-selected cross-validation points and their indicies. 
#' 
#' @param modelType Method to be used for model-fitting. Valid options are: \code{glmnet},\code{pamr}, and \code{sam}.
#' @param thresholdCVRates Matrix of error rates at regularizationThresholds returned by \code{citrus.thresholdCVs.*} function.
#' @param fdrRate FDR Maximum used to determine FDR-constrained model regularization threshold.
#' 
#' @return List of regularization thresholds and indicies based on pre-selected cross-validation error rates points.
#' 
#' @details For predictive models (i.e. \code{pamr} or \code{glmnet}), returns indicies of regularization thresholds
#' producing the minimum cross validation error rate (\code{cv.min}), the simplest model having error within 1
#' standard error of the minimum (\code{cv.1se}), and the model with the minimum error having an FDR rate < \code{fdrRate} (\code{cv.fdr.constrained})
#' when possible.
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Calculate regularization thresholds
#' regularizationThresholds = citrus.generateRegularizationThresholds(abundanceFeatures,labels,modelType="pamr",family="classification")
#' 
#' # Calculate CV Error rates
#' thresholdCVRates = citrus.thresholdCVs.quick("pamr",abundanceFeatures,labels,regularizationThresholds,family="classification") 
#' 
#' # Get pre-selected CV Minima
#' cvMinima = citrus.getCVMinima("pamr",thresholdCVRates)
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

#' Report model features at pre-specified thresholds.
#' 
#' Report model features at pre-specific thresholds. For predictive models, reports non-zero model features at specified
#' regularization thresholds. For FDR-constrained models, reports features below specified false discovery rates.
#' 
#' @param cvMinima List of regularization indicies at which to extract model features, produced by \code{\link{citrus.getCVMinima}}.
#' @param finalModel Predictive model from which to extract non-zero features.
#' @param finalFeatures Features used to construct \code{finalModel}.
#' 
#' @return List of significant features and clusters at specified thresholds.
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' abundanceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,clusterAssignments=citrus.clustering$clusterMembership,clusterIds=largeEnoughClusters)
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Calculate regularization thresholds
#' regularizationThresholds = citrus.generateRegularizationThresholds(abundanceFeatures,labels,modelType="pamr",family="classification")
#' 
#' # Calculate CV Error rates
#' thresholdCVRates = citrus.thresholdCVs.quick("pamr",abundanceFeatures,labels,regularizationThresholds,family="classification") 
#' 
#' # Get pre-selected CV Minima
#' cvMinima = citrus.getCVMinima("pamr",thresholdCVRates)
#' 
#' # Build Final Model
#' finalModel = citrus.buildEndpointModel(abundanceFeatures,labels,family="classification",type="pamr",regularizationThresholds)
#' 
#' # Get model features
#' citrus.extractModelFeatures(cvMinima,finalModel,abundanceFeatures)
citrus.extractModelFeatures = function(cvMinima,finalModel,finalFeatures){
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

calculatePredictionErrorRate = function(predictionScore,regularizationThresholds,family){
  nFolds=length(predictionScore)
  counter=1;
  tmp=list()
  for (i in 1:nFolds){
    for (j in 1:nrow(predictionScore[[i]])){
      tmp[[counter]] = predictionScore[[i]][j,]
      length(tmp[[counter]])=length(regularizationThresholds)
      counter=counter+1;
    }
  }
  bound = do.call("rbind",tmp)
  thresholdMeans = apply(bound,2,mean,na.rm=T)
  
  thresholdSEMs = apply(bound,2,sd,na.rm=T)/sqrt(apply(!is.na(bound),2,sum))
  return(list(cvm=thresholdMeans,cvsd=thresholdSEMs))
} 
