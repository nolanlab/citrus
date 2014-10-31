#' Perform a full Citrus analysis
#' 
#' This function takes FCS data, identifies cell subsets using clustering, 
#' characterizes the behavior of each identified cluster on a per-sample basis,
#' and identifies cluster properties associated with an endpoint.
#' 
#' @param fileList A data frame containing the names of FCS files to be analyzed. Each column should contain FCS files measured in the same condition. Each row should contain FCS files measured from the same individual under different conditions. See examples below.
#' @param labels A vector containing the experimental endpoint values for each row of FCS files in the fileList.
#' @param clusteringColumns A vector containing the names or indices of parameters to be used for clustering.
#' @param dataDirectory The full path to the directory containing FCS files.
#' @param outputDirectory The full path to directory analysis output should be placed in. If not \code{NULL}, clustering is saved.
#' @param family Type of association model to be calculated. Valid values are currently \code{classification} and \code{continuous}.
#' @param modelTypes Vector of model types to be used to detect associations with experimental endpoint. Valid values are \code{glmnet}, \code{pamr}, and \code{sam}.
#' @param nFolds Number of independent clustering folds to be used for model fitting. Should only be >1 if using \code{glmnet} or \code{pamr} models. Default value is 1 and all samples are clustered together.
#' @param conditionComparaMatrix matrix of condition data to compare. See details.
#' @param ... Other arguments to be passed to \code{citrus.readFCSSet}, \code{citrus.calculateFeatures}, \code{citrus.endpointRegress}. 
#' Useful additional arguments listed in details.
#' 
#' @details Citrus is able to analyze FCS data from single conditions or relative to a given baseline 
#' (i.e. stimulated values relative to an unstimulated baseline). Tell Citrus to compare conditions using 
#' the conditionComparisonMatrix argument. The conditionComparisonMatrix is an \eqn{n x n} matrix of logical values
#' indicating which conditions should be compared where \eqn{n} is the number of experimental conditions specified in
#' the fileList. Rows and columns names of the conditionComparsion matrix should be the names of experimental conditions 
#' in the fileList. A matrix entry of \code{TRUE} indicates that Citrus should calculate the value features in column name
#' relative to feature values in the row name. \code{TRUE} entries on the diagonal means conditions are analyzed by themselves. 
#' See Examples.
#' 
#' Other useful arguments to pass to citrus.full include:
#' \itemize{
#' 
#' \item \code{fileSampleSize}: The number of cells to select from each sample. See \code{\link{citrus.readFCSSet}}.
#' 
#' \item \code{transformColumns}: Vector of parameter names or indicies whose values should transformed. See \code{\link{citrus.readFCSSet}}.
#' 
#' \item \code{transformCofactor}:  Cofactor for arcsin-hyperbolic transform. See \code{\link{citrus.readFCSSet}}.
#' 
#' \item \code{featureType}: The descriptive feature type to be calculated for each cluster. See \code{\link{citrus.calculateFeatures}}.
#' 
#' \item \code{minimumClusterSizePercent}: The minimum cluster size as a percentage of total sampled cells. See \code{\link{citrus.calculateFeatures}}.
#' }
#' 
#' @return A citrus.full.result object with the following properties:
#' \item{citrus.combinedFCSSet}{A \code{citrus.combinedFCSSet} object containing data read from files.}
#' \item{citrus.foldClustering}{A \code{citrus.foldClustering} object that contains clustering results.}
#' \item{conditionFoldFeatures}{A list of \code{citrus.foldFeatureSet} objects used for regression, one list element for each set of analyzed conditions.}
#' \item{conditionRegressionResults}{A list of \code{citrus.regressionResult} objects, one entry for each set of analyzed condition.}
#' \item{family}{The family of the regression model.}
#' \item{labels}{Endpoint labels of the analyzed samples.}
#' \item{conditions}{List of analyzed conditions.}
#' \item{modelTypes}{Vector of regression models used to test for endpoint associations.}
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @seealso \code{citrus.readFCSSet}, \code{citrus.cluster},\code{citrus.calculateFeatures}, and \code{citrus.endpointRegress}
#' 
#' @examples
#' # Example with a single experimental condition (unstimulated data) 
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Run citrusAnalysis 
#' results = citrus.full(fileList,labels,clusteringColumns,dataDirectory)
#' 
#' # Should be used to plot results
#' # plot(results, outputDirectory="/path/to/outputDirectory")
#' 
#' ############################################
#' 
#' # Example using multiple experimental conditions, including a 
#' # comparison of stimulated to unstimulated damples
#' 
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example2")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame(unstim=list.files(dataDirectory,pattern="unstim"),stim1=list.files(dataDirectory,pattern="stim1"))
#' 
#' # Disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Vector of parameters to be used for clustering
#' clusteringColumns = c("LineageMarker1","LineageMarker2")
#' 
#' # Vector of parameters to calculate medians for
#' functionalColumns = c("FunctionalMarker1","FunctionalMarker2")
#' 
#' # Form condition comparison matrix to tell Citrus which conditions to analyze compare
#' conditionComparaMatrix = matrix(T,nrow=2,ncol=2,dimnames=list(colnames(fileList),colnames(fileList)))
#' conditionComparaMatrix[2]=F
#' 
#' # Run citrusAnalysis 
#' results = citrus.full(fileList,labels,clusteringColumns,dataDirectory,
#'                       conditionComparaMatrix=conditionComparaMatrix,
#'                       featureType="medians",
#'                       medianColumns=functionalColumns)
#'
#' # Result should be plotted
#' # plot(results,outputDirectory="/path/to/outputDirectory")
citrus.full = function(fileList,
                       labels,
                       clusteringColumns,
                       dataDirectory,
                       outputDirectory=NULL,
                       family="classification",
                       modelTypes=c("glmnet"),
                       nFolds=1,
                       conditionComparaMatrix=NULL,
                       ...){
  
  # No point in running cv if SAM only model
  if (all(modelTypes=="sam")){
    nFolds=1
  }
  
  if (!is.null(conditionComparaMatrix)){
    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix)
  } else {
    allConditions = as.list(colnames(fileList))
  }
  names(allConditions) = sapply(allConditions,function(x){paste(rev(x),collapse="_vs_")})
  
  # Read in data
  citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory=dataDirectory,fileList=fileList,...)
    
  # Cluster each fold
  citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,labels=labels,nFolds=nFolds)
  if (!is.null(outputDirectory)){
    save(list=c("citrus.combinedFCSSet","citrus.foldClustering"),file=file.path(outputDirectory,"citrusClustering.rData"),compress=F)  
  }
  
    
  conditionFeatures = list()
  conditionRegressionResults = list()
  # Analyze each condition or set of conditions separately 
  for (conditions in allConditions){
    cat(paste0("Analyzing condition(s): ",paste0(conditions,collapse=" vs. "),"\n"))
    #conditionFileIds = citrus.combinedFCSSet$fileIds[,conditions]
    
    # Calculate fold features
    #citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,minimumClusterSizePercent=minimumClusterSizePercent,medianColumns=medianColumns,mc.cores=4)
    #citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,conditions=conditions,featureType="medians",medianColumns=medianColumns)
    cat("\tBuilding Fold Features\n")
    citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering=citrus.foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,conditions=conditions,...)
    conditionFeatures[[paste(rev(conditions),collapse="_vs_")]] = citrus.foldFeatureSet
    
    # Endpoint regress for each model type
    #citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
    #citrus.regressionResults = citrus.endpointRegress("sam",citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family)
    cat("\tAnalyzing vs. endpoint\n")
    citrus.regressionResults = mclapply(modelTypes,citrus.endpointRegress,citrus.foldFeatureSet=citrus.foldFeatureSet,labels=labels,family=family,...)
    names(citrus.regressionResults) = modelTypes
    conditionRegressionResults[[paste(rev(conditions),collapse="_vs_")]] = citrus.regressionResults
    cat("\n")
  }
  # Return Results
  results = list(citrus.combinedFCSSet=citrus.combinedFCSSet,citrus.foldClustering=citrus.foldClustering,conditionFoldFeatures=conditionFeatures,conditionRegressionResults=conditionRegressionResults,family=family,labels=labels,conditions=allConditions,modelTypes=modelTypes)
  class(results) = "citrus.full.result"
  return(results)
}

#' @export
#' @name citrus.full
print.citrus.full.result = function(citrus.full.result,...){
  cat("Citrus full result\n")
  cat(paste("Family:",citrus.full.result$family,"\n"))
  cat("Models: ")
  cat(paste0(paste(unlist(citrus.full.result$modelTypes),collapse=", "),"\n"))
  cat("Conditions: ")
  cat(paste0(paste(unlist(citrus.full.result$conditions),collapse=", "),"\n"))
}

