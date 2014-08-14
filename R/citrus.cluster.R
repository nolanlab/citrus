#' Map new data to existing clusters
#' 
#' Map new data to clusters defined by a clustering. 
#' 
#' @param citrus.combinedFCSSet.new A \code{citrus.combinedFCSSet} object containing new data to be mapped.
#' @param citrus.combinedFCSSet.old A \code{citrus.combinedFCSSet} object containing new data that has been clustered.
#' @param citrus.clustering A clustering of data in \code{citrus.combinedFCSSet.old}.
#' @param mappingColumns Parameters to be used for mapping new data to existing cluster space. If \code{NULL}, then parameters 
#' used to cluster \code{citrus.combinedFCSSet.old} are used.
#' @param ... Additional arguments (ignored).
#' 
#' @return A \code{citrus.mapping} object
#' \item{clusterMembership}{List of event indices, reporting which events in \code{citrus.combinedFCSSet.new} belong to each cluster in \code{citrus.clustering}.}
#' \item{mappingColumns}{Columns used to match new data events with existing clustered data.}
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
#' citrus.mapping = citrus.mapToClusterSpace(citrus.combinedFCSSet.new=citrus.combinedFCSSet1,citrus.combinedFCSSet.old=citrus.combinedFCSSet2,citrus.clustering)
citrus.mapToClusterSpace = function(citrus.combinedFCSSet.new,citrus.combinedFCSSet.old,citrus.clustering,mappingColumns=NULL,...){
  if (is.null(mappingColumns)){
    mappingColumns = citrus.clustering$clusteringColumns
  }
  
  # Map new data to nearest neighbor in existing clustered dataset
  cat(paste("Mapping",nrow(citrus.combinedFCSSet.new$data),"events\n"))
  nearestNeighborMap = citrus.assignToCluster(tbl=citrus.combinedFCSSet.new$data[,mappingColumns],cluster_data=citrus.combinedFCSSet.old$data[,mappingColumns],cluster_assign=rep(1,nrow(citrus.combinedFCSSet.old$data)))
  
  # Assign new data to the same clusters as its nearest neighbor in the existing clustered dataset
  newDataClusterAssignments = mclapply(citrus.clustering$clusterMembership,citrus.mapNeighborsToCluster,nearestNeighborMap=nearestNeighborMap,...)
  
  result = list(clusterMembership=newDataClusterAssignments,mappingColumns=mappingColumns)
  class(result) = "citrus.mapping"
  return(result)
}

citrus.mapNeighborsToCluster = function(clusterMembers,nearestNeighborMap){
  which(nearestNeighborMap %in% clusterMembers)
}

citrus.mapFoldDataToClusterSpace = function(foldIndex,folds,foldClustering,citrus.combinedFCSSet,...){
  leavoutFileIds = as.vector(citrus.combinedFCSSet$fileIds[folds[[foldIndex]],])
  includeFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),leavoutFileIds)
  
  
  citrus.mapToClusterSpace(citrus.combinedFCSSet.new=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,leavoutFileIds),
                           citrus.combinedFCSSet.old=citrus.maskCombinedFCSSet(citrus.combinedFCSSet,includeFileIds),
                           citrus.clustering=foldClustering[[foldIndex]],
                           mappingColumns=foldClustering[[foldIndex]]$clusteringColumns,...)
}


citrus.traverseMergeOrder = function(node,mergeOrder){
  left = mergeOrder[node,1]
  right = mergeOrder[node,2]
  if (left<0){
    leftAtom = -left;  
  } else {
    leftAtom = citrus.traverseMergeOrder(left,mergeOrder);
  }
  
  if (right<0){
    rightAtom = -right;
  } else {     
    rightAtom = citrus.traverseMergeOrder(right,mergeOrder); 
  }
  return(c(leftAtom,rightAtom));
}

#' Get ancestors or decendants of a cluster
#'
#' Get ancestors or decendants of a cluster for hierarchical clustering
#' @name citrus.getRelatedClusterIds
#' @param clusterId ID of a cluster for which to retreive related clusters. 
#' @param mergeOrder \code{mergeOrder} result from \code{hclust} object.
#' 
#' @return Vector of cluster ids that are decendants of \code{clusterId} argument.
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
#' # Get decendants
#' citrus.getClusterDecendants(15000,citrus.clustering$clustering$merge)
#' 
#' # Get ancestors
#' citrus.getClusterAncestors(15000,citrus.clustering$clustering$merge)
citrus.getClusterDecendants = function(clusterId,mergeOrder){
  left = mergeOrder[clusterId,1]
  right = mergeOrder[clusterId,2]
  if (left>0){
    left = c(left,citrus.getClusterDecendants(left,mergeOrder))
  } else {
    left = c()
  }
  if (right>0){
    right = c(right,citrus.getClusterDecendants(right,mergeOrder))
  } else {
    right = c()
  }
  return(c(left,right))
}

#' @rdname citrus.getRelatedClusterIds
#' @name citrus.getRelatedClusterIds
#' @export
citrus.getClusterAncestors = function(clusterId,mergeOrder){
  parent = which(mergeOrder==clusterId,arr.ind=T)[1]
  if (is.na(parent)){
    return(c())
  } else {
    return(c(parent,citrus.getClusterAncestors(parent,mergeOrder)))
  }
}

#' Cluster a \code{citrus.combinedFCSSet}
#' 
#' Cluster data in a \code{citrus.combinedFCS} set object
#' 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param clusteringColumns Vector of names or indicies of data to be used for clustering.
#' @param clusteringType Type of clustering to be perfomed. Valid options are: \code{hierarchical}.
#' @param ... Other arguments passed to specific clustering algorithm.
#' 
#' @return A \code{citrus.clustering} object.
#' \item{clustering}{Clustering object returned by clustering algorithm.}
#' \item{clusterMembership}{List, containing indicies of cells in \code{citrus.combinedFCSSet} data that belong to each cluster.}
#' \item{type}{Type of clustering performed.}
#' \item{clusteringColumns}{Parameters used for clustering}
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
citrus.cluster = function(citrus.combinedFCSSet,clusteringColumns,clusteringType="hierarchical",...){
  clustering = do.call(paste0("citrus.cluster.",clusteringType),args=list(data=citrus.combinedFCSSet$data[,clusteringColumns]))
  clusterMembership = do.call(paste0("citrus.calculateClusteringMembership.",clusteringType),args=list(clustering=clustering))
  result = list(clustering=clustering,clusterMembership=clusterMembership,type=clusteringType,clusteringColumns=clusteringColumns)
  class(result) = "citrus.clustering"
  return(result)
}


#' @export
#' @name citrus.cluster
print.citrus.clustering = function(citrus.clustering){
  cat("Citrus clustering\n")
  cat(paste0("\tType:\t\t",citrus.clustering$type,"\n"))
  cat(paste0("\tClusters:\t",length(citrus.clustering$clusterMembership),"\n"))
}

citrus.clusterFold = function(foldIndex,folds,citrus.combinedFCSSet,clusteringColumns,...){
  leavoutFileIds = as.vector(citrus.combinedFCSSet$fileIds[folds[[foldIndex]],])
  includeFileIds = setdiff(as.vector(citrus.combinedFCSSet$fileIds),leavoutFileIds)
  citrus.cluster(citrus.maskCombinedFCSSet(citrus.combinedFCSSet,fileIds=includeFileIds),clusteringColumns,...)
}

#' Cluster independent folds of data
#' 
#' Cluster subsets of data from different samples and maps leftout sample data to fold cluster space.
#' 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param clusteringColumns Vector of parameter names or indicies to be used for clustering.
#' @param labels Labels of samples being clustered. If supplied, used for balancing folds for clustering
#' @param nFolds Number of independent folds of clustering to perform. If \code{nFolds=1}, all data are 
#' clustered together and model is regression model is constructed from single feature set.
#' @param ... Other arguments passed to specific clustering functions.
#' 
#' @return A \code{citrus.foldClustering} object
#' \item{folds}{Indicies of sample rows to be omitted during each fold of clustering. Only defined if \code{nFolds > 1}.}
#' \item{foldClustering}{\code{citrus.clustering} objects for each fold. Only defined if \code{nFolds > 1}.}
#' \item{foldMappingAssignments}{\code{citrus.mapping} objects for each fold, containing mapping of data from left-out samples. Only defined if \code{nFolds > 1}}
#' \item{allClustering}{\code{citrus.clusteringObject} from all sample data.}
#' \item{nFolds}{Number of independent folds.}
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
citrus.clusterAndMapFolds = function(citrus.combinedFCSSet,clusteringColumns,labels=NULL,nFolds=1,...){
  
  result = list()
  
  if (nFolds>1){
    
    # This should eventually changed to just make fold w/o balancing.
    if (is.null(labels)){
      stop("Labels must be supplied for balanced fold selection.")
    }
    
    # Define Folds
    if (is.numeric(labels)){
      result$folds = pamr:::balanced.folds(nfolds=nFolds,y=sample(rep(1:nFolds,length.out=length(labels))))
    } else {
      result$folds = pamr:::balanced.folds(y=labels,nfolds=nFolds)  
    }
    
    
    # FOLDS MUST BE SORTED
    result$folds = lapply(result$folds,sort)
    
    # Cluster each fold of data together
    result$foldClustering = lapply(1:nFolds,citrus.clusterFold,folds=result$folds,citrus.combinedFCSSet=citrus.combinedFCSSet,clusteringColumns=clusteringColumns,...)
    
    # Map the left out data to the fold clustered data
    result$foldMappingAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,folds=result$folds,foldClustering=result$foldClustering,citrus.combinedFCSSet=citrus.combinedFCSSet)  
  }
  
  # Cluster all the data together
  result$allClustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns,...)
  
  result$nFolds = nFolds
  
  class(result) = "citrus.foldClustering"
  return(result)
}

#' @export
#' @name citrus.clusterAndMapFolds
print.citrus.foldClustering = function(citrus.foldClustering){
  cat("citrus.foldClustering\n")
  cat(paste("\tNumber of Folds: ",length(citrus.foldClustering$folds),"\n"))
}

citrus.cluster.hierarchical = function(data){
  cat(paste("Clustering",nrow(data),"events\n"));
  return(Rclusterpp.hclust(data))  
}

citrus.calculateClusteringMembership.hierarchical = function(clustering,...){
  clusterType="hierarchical"
  if ("clusterType" %in% names(list(...))){
    clusterType = list(...)[["clusterType"]]
  }
  if (clusterType=="hierarchical"){
    mergeOrder = clustering$merge
    if ("mc.cores" %in% names(list(...))){
      return(mclapply(as.list(1:nrow(mergeOrder)),citrus.traverseMergeOrder,mergeOrder=mergeOrder,mc.cores=list(...)[["mc.cores"]]))
    } else {
      return(lapply(as.list(1:nrow(mergeOrder)),citrus.traverseMergeOrder,mergeOrder=mergeOrder))
    }
  } else {
    stop(paste("Don't know how to caclulate complete clustering for clusterType",clusterType))
  }
}