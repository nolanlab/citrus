#' Calculate descriptive cluster features
#' 
#' Calculate descriptive properties for each cluster in each sample
#' 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param clusterAssignments List of indicies of cluster assignments for each cluster.
#' @param clusterIds Vector of cluster ID's for which to calculate features.
#' @param featureType Type of feature to calculate. Valid options are: \code{abundances} and \code{medians}.
#' See \code{citrus.calculateFeature.*} functions for feature-type specific arguments. 
#' @param conditions Vector of conditions for which to calculate features. See details.
#' @param ... Other arguments passed to individual feature-calculation functions.
#' 
#' @details If \code{conditions=NULL}, \code{citrus.calculateFeatures} constructs features for all samples 
#' in the \code{citrus.combinedFCSSet}. If \code{conditions} is a single element, \code{citrus.calculateFeatures} 
#' constructs features for samples in that condition. If \code{conditions} contains two elements, the first
#' condition is used as a baseline condition, the second is used as a comparison condition, and \code{citrus.calculateFeatures}
#' returns the difference in feature values between the comparison and baseline conditions. 
#' 
#' @return Matrix of cluster features
#' 
#' @seealso \code{citrus.calculateFeature.type}
#' @author Robert Bruggner
#' @export 
#' 
#' @examples
#' ######################################################
#' # Calculate cluster abundances for single condition
#' ######################################################
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
#' 
#' ######################################################
#' # Calculate median levels of functional markers in 
#' # stimulated conditions relative to unstimluated
#' # condtion 
#' ######################################################
#' # Where the data Lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example2")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame(unstim=list.files(dataDirectory,pattern="unstim"),stim1=list.files(dataDirectory,pattern="stim1"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # Vector of parameters to be used for clustering
#' clusteringColumns = c("LineageMarker1","LineageMarker2")
#' 
#' # Vector of parameters to calculate medians for
#' functionalColumns = c("FunctionalMarker1","FunctionalMarker2")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Build features
#' medianDifferenceFeatures = citrus.calculateFeatures(citrus.combinedFCSSet,
#'                                                 clusterAssignments=citrus.clustering$clusterMembership,
#'                                                 clusterIds=largeEnoughClusters,
#'                                                 featureType="medians",
#'                                                 medianColumns=functionalColumns,
#'                                                 conditions=c("unstim","stim1"))
citrus.calculateFeatures = function(citrus.combinedFCSSet,clusterAssignments,clusterIds,featureType="abundances",conditions=NULL,...){  

  fileIds=NULL
  baselineFileIds=NULL
  
  # If conditions are provided, only calculate features for condition files.
  if (!is.null(conditions)){
    # If one condition, just calculate features for that condition.
    if (length(conditions)==1){
      fileIds = citrus.combinedFCSSet$fileIds[,conditions[1]]
      # If two conditions are provided, calculate the difference between those conditions
    } else if (length(conditions)==2){
      fileIds = citrus.combinedFCSSet$fileIds[,conditions[2]]
      baselineFileIds = citrus.combinedFCSSet$fileIds[,conditions[1]]
    } else {
      stop("Unexpected number of conditions")
    }
  }
  
  #Error check before we actually start the work.
  if ((!all(featureType %in% citrus.featureTypes()))||(length(featureType)<1)){
    stop(paste("featureType must be one of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
  }
  
  # If provided, calculate the value of a feature relative to some baseline.
  if ((!is.null(fileIds)) && (!is.null(baselineFileIds))){
    baselineFeatures = do.call(paste0("citrus.calculateFeature.",featureType),args=list(clusterIds=clusterIds,clusterAssignments=clusterAssignments,citrus.combinedFCSSet=citrus.combinedFCSSet,fileIds=baselineFileIds,...))  
    referenceFeatures = do.call(paste0("citrus.calculateFeature.",featureType),args=list(clusterIds=clusterIds,clusterAssignments=clusterAssignments,citrus.combinedFCSSet=citrus.combinedFCSSet,fileIds=fileIds,...))  
    differenceFeatures = referenceFeatures-baselineFeatures
    colnames(differenceFeatures) = paste(colnames(differenceFeatures),"difference")
    return(differenceFeatures)
  } else {
    # Otherwise, just calculate the features for every available file.  
    return(do.call(paste0("citrus.calculateFeature.",featureType),args=list(clusterIds=clusterIds,clusterAssignments=clusterAssignments,citrus.combinedFCSSet=citrus.combinedFCSSet,fileIds=fileIds,...))  )
  }
  
}

#' Calculate descriptive cluster features
#' 
#' Calculate descriptive cluster features.
#' 
#' @name citrus.calculateFeature.type
#' @param clusterIds Cluster IDs that descriptive features should be calculated for.
#' @param clusterAssignments List with indicies of cells belonging to each cluster. 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param fileIds Vector of file IDs to calculate features for. If \code{NULL}, calculates
#' features for all samples in \code{citrus.combinedFCSSet}.
#' @param medianColumns Vector of parameter names or numeric indicies of parameters for which to calculate cluster median values for. 
#' @param ... Other arguments (ignored).
#' 
#' @details See \code{citrus.calculateFeatures} for examples.
#' 
#' @author Robert Bruggner
#' @export 
#' 
#' @seealso \code{citrus.calculateFeatures}, \code{citrus.calculateFoldFeatureSet}
citrus.calculateFeature.abundances = function(clusterIds,clusterAssignments,citrus.combinedFCSSet,fileIds=NULL,...){
  eventFileIds = citrus.combinedFCSSet$data[,"fileId"]
  if (is.null(fileIds)){
    fileIds = unique(eventFileIds)
  }
  res = mcmapply(citrus.calculateFileClusterAbundance,
          clusterId=rep(clusterIds,length(fileIds)),
          fileId=rep(fileIds,each=length(clusterIds)),
          MoreArgs=list(
                        clusterAssignments=clusterAssignments,
                        eventFileIds=eventFileIds))
  return(
    matrix(res,ncol=length(clusterIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames[fileIds],paste("cluster",clusterIds,"abundance")))
  )
}


citrus.calculateFileClusterAbundance = function(clusterId,fileId,clusterAssignments,eventFileIds,...){
  res = sum(which(eventFileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((eventFileIds==fileId))
  # Multiply by 100 so on a scale from 0-100 instead of 0-1.
  res*100
}

# Median Features
#' @rdname citrus.calculateFeature.type
#' @name citrus.calculateFeature.type
#' @export
citrus.calculateFeature.medians = function(clusterIds,clusterAssignments,citrus.combinedFCSSet,medianColumns,fileIds=NULL,...){
  
  if (is.null(fileIds)){
    fileIds = unique(citrus.combinedFCSSet$data[,"fileId"])  
  }
  
  res = mcmapply(citrus.calculateFileClusterMedian,
                 clusterId=rep(rep(clusterIds,length(medianColumns)),length(fileIds)),
                 fileId=rep(fileIds,each=(length(clusterIds)*length(medianColumns))),
                 medianColumn=rep(rep(medianColumns,each=length(clusterIds)),length(fileIds)),
                 MoreArgs=list(
                   data=citrus.combinedFCSSet$data,
                   clusterAssignments=clusterAssignments))
  
  return(
    matrix(res,nrow=length(fileIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames[fileIds],paste("cluster",rep(clusterIds,length(medianColumns)),rep(medianColumns,each=length(clusterIds)),"median")))
  )
}

citrus.calculateFileClusterMedian = function(clusterId,fileId,medianColumn,data,clusterAssignments,...){
  if (length(clusterAssignments[[clusterId]])<3){
    return(0)
  }
  clusterData = data[clusterAssignments[[clusterId]],]
  if (sum(clusterData[,"fileId"]==fileId)<3){
    return(0)
  }
  clusterFileDataValues = clusterData[clusterData[,"fileId"]==fileId,medianColumn]
  return(median(clusterFileDataValues))
}


###############################################
# Functions for dealing with fold features
###############################################
#' Build cluster features for folds of clustering
#'
#' Build cluster features for each fold of clustering. If multiple folds of clustering have been performed, \code{citrus.calculateFoldFeatureSet}
#' builds features for clustered and leftout samples for each fold.
#' 
#' @param citrus.foldClustering A \code{citrus.foldClustering} object
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object
#' @param featureType Type of feature to be calculated. Valid options are: \code{abundances} and \code{medians}. See \code{\link{citrus.calculateFeatures}} for additional argument details.
#' @param minimumClusterSizePercent Minimum cluster size percent used to select clusters for analysis. See \code{\link{citrus.selectClusters}}.
#' @param ... Additional arguments passed to feature-type specific calculation functions.
#' 
#' @return A \code{citrus.foldFeatureSet} object with properties:
#' \item{foldLargeEnoughClusters}{List of selected clusters for each fold of clustering.}
#' \item{foldFeatures}{List of features constructed from fold clustered samples.}
#' \item{leftoutFeatures}{List of features constructed from non-clustered samples that were mapped to the fold clustering space.}
#' \item{allLargeEnoughClusters}{Selected clusters from clustering of all samples.}
#' \item{allFeatures}{Features constructed from clustering of all samples.}
#' \item{minimumClusterSizePercent}{User-specified minimum cluster size percent.}
#' \item{folds}{List of sample folds.}
#' \item{nFolds}{Number of folds.}
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
citrus.calculateFoldFeatureSet = function(citrus.foldClustering,citrus.combinedFCSSet,featureType="abundances",minimumClusterSizePercent=0.05,...){
  result = list()
    
  if (citrus.foldClustering$nFolds>1){
    # Select clusters in each folds
    # Default is minimum cluster size
    result$foldLargeEnoughClusters = lapply(citrus.foldClustering$foldClustering,citrus.selectClusters,minimumClusterSizePercent=minimumClusterSizePercent)  
    
    # Build Training Features
    #result$foldFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.calculateFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType="abundances",foldClustering=citrus.foldClustering$foldClustering,conditions=conditions)
    result$foldFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.calculateFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,foldClustering=citrus.foldClustering$foldClustering,...)
    
    # Build Testing Features
    #result$leftoutFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.calculateFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType="abundances",foldMappingAssignments=citrus.foldClustering$foldMappingAssignments,calculateLeftoutFeatureValues=T,conditions=conditions)
    result$leftoutFeatures = lapply(1:citrus.foldClustering$nFolds,citrus.calculateFoldFeatures,folds=citrus.foldClustering$folds,foldClusterIds=result$foldLargeEnoughClusters,citrus.combinedFCSSet=citrus.combinedFCSSet,featureType=featureType,foldMappingAssignments=citrus.foldClustering$foldMappingAssignments,calculateLeftoutFeatureValues=T,...)
  }
  
  # Build features for clustering of all samples
  result$allLargeEnoughClusters = citrus.selectClusters(citrus.clustering=citrus.foldClustering$allClustering,minimumClusterSizePercent=minimumClusterSizePercent)
  #result$allFeatures = citrus.calculateFeatures(citrus.combinedFCSSet=citrus.combinedFCSSet,clusterAssignments=citrus.foldClustering$allClustering$clusterMembership,clusterIds=result$allLargeEnoughClusters,featureType="abundances",conditions=conditions)
  result$allFeatures = citrus.calculateFeatures(citrus.combinedFCSSet=citrus.combinedFCSSet,clusterAssignments=citrus.foldClustering$allClustering$clusterMembership,clusterIds=result$allLargeEnoughClusters,featureType=featureType,...)
  
  # Extra feature building parameters, etc
  result$minimumClusterSizePercent=minimumClusterSizePercent
  result$folds=citrus.foldClustering$folds
  result$nFolds=citrus.foldClustering$nFolds
  
  class(result) = "citrus.foldFeatureSet"
  return(result)
}

citrus.calculateFoldFeatures = function(foldIndex,folds,foldClusterIds,citrus.combinedFCSSet,foldClustering=NULL,foldMappingAssignments=NULL,featureType="abundances",calculateLeftoutFeatureValues=F,...){
  
  leavoutFileIndices = folds[[foldIndex]]
    
  if (calculateLeftoutFeatureValues){
    featureFileIds = as.vector(citrus.combinedFCSSet$fileIds[leavoutFileIndices,])
    if (is.null(foldMappingAssignments)){
      stop("foldMappingAssignments argument missing")
    }
    cellAssignments = foldMappingAssignments[[foldIndex]]$clusterMembership
  } else {
    featureFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),as.vector(citrus.combinedFCSSet$fileIds[leavoutFileIndices,]))
    if (is.null(foldClustering)){
      stop("foldClustering argument missing")
    }
    cellAssignments = foldClustering[[foldIndex]]$clusterMembership
  }
  
  citrus.calculateFeatures(citrus.combinedFCSSet=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,featureFileIds),
                       clusterAssignments=cellAssignments,
                       clusterIds=foldClusterIds[[foldIndex]],
                       featureType=featureType,
                       ...)
  
}