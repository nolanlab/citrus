#' Read a set of FCS files for analysis
#' 
#' Reads and combines data from many FCS files for use with Citrus. 
#' 
#' @param dataDirectory Full or relative path to directory containin FCS files.
#' @param fileList Data frame containing names of files to be read. Each entry should be a file name. FCS files 
#' from the same experimental condition should be listed in the same column with the column name indicating the experimental 
#' condition. FCS files from the same sample measured in multiple experimental conditions should be listed in the same row.
#' See examples. 
#' @param fileSampleSize Number of cells to be selected from each FCS file. If this number is larger than the number of events
#' in a given FCS file, all cells from that file are selected. 
#' @param transformColumns Vector of parameter names or indicies whose values should transformed prior to analysis.
#' @param Cofactor for arcsin-hyperbolic transform.
#' @param scaleColumns Vector of parameter names or indicies whose values should be scaled prior to analysis.
#' @param useChannelDescriptions Should channel descriptive names be used to name data instead of channel names?
#' @param readParameters Vector of parameter names (or descriptions if \code{useChannelDescriptions='T'}) to be read from each FCS file. Can be used to specify common parameters
#' among FCS files that contain different measured parameters.
#' @param ... Other undocumented parameters.
#' 
#' @return An object of type \code{citrus.combinedFCSSet}.
#' \item{data}{Combined FCS measurements from files listed in fileList, with paraemters \code{fileEventNumber} and \code{fileId} added. The former
#' reports the event number in the original source file and the latter records the file from which the event was sampled.}
#' \item{fileIds}{A matrix recording the assigned fileIds for each file in the fileList argument.}
#' \item{fileNames}{The corresponding file names for each file, indexed by fileId.}
#' \item{fileChannelNames}{Names of parameters in each read file.}
#' \item{fileReagentNames}{Description of parameters in each read file.}
#' \item{transformColumns}{List of columns that were transformed when read, if supplied.}
#' \item{transformCofactor}{Cofactor used for \code{arcsinh} transformation, if applicable.}
#' 
#' @author Robert Bruggner
#' @export
#'  
#' @examples
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed - one condition.
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # Create list of files to be analyzed - two conditions.
#' fileList = data.frame(unstim=list.files(dataDirectory,pattern="unstim"),stim1=list.files(dataDirectory,pattern="stim1"))
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
citrus.readFCSSet = function(dataDirectory,fileList,fileSampleSize=1000,transformColumns=NULL,transformCofactor=5,scaleColumns=NULL,useChannelDescriptions=F,readParameters=NULL,...){
  data = list();
  fileCounter = 1;
  fileNames = c();
  fileChannelNames = list();
  fileReagentNames = list();
  addtlArgs = list(...)
    
  conditions = colnames(fileList)
  
  for (i in 1:length(conditions)){
    cat(paste("Reading Condition ",conditions[i],"\n"));
    conditionData = list();
    fileChannelNames[[conditions[i]]] = list();
    fileReagentNames[[conditions[i]]] = list();
    for (fileName in fileList[,conditions[i]]){
      fileNames[fileCounter]=fileName
      filePath = file.path(dataDirectory,fileName);
      if (!file.exists(filePath)){
        stop(paste("File",filePath,"not found."));
      }
      cat(paste("\tReading file ",fileName,"\n"));
      fcsFile = citrus.readFCS(filePath)
      
      fcsData=exprs(fcsFile)
      
      if (useChannelDescriptions){
       channelDescriptions = as.vector(pData(parameters(fcsFile))$desc)
       colnames(fcsData)[nchar(channelDescriptions)>2] = channelDescriptions[nchar(channelDescriptions)>2]
      }
      
      if (!is.null(readParameters)){
        fcsData = fcsData[,readParameters]
      }
      
      fileChannelNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$name)
      fileReagentNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$desc)
      fcsData = cbind(fcsData,fileEventNumber=1:nrow(fcsData),fileId=fileCounter);
      fileCounter=fileCounter+1;
      
      if ((!is.null(fileSampleSize))&&(fileSampleSize<nrow(fcsData))){
        cat(paste("\tSampling",fileSampleSize,"events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),fileSampleSize)),] 
      }
      
      conditionData[[fileName]] = fcsData
      
      if (useChannelDescriptions){
        channelDescriptions = as.vector(pData(parameters(fcsFile))$desc)
        nchar(channelDescriptions)>2
      }
    }
    data[[conditions[i]]] = do.call("rbind",conditionData)
    rm(conditionData)
    gc();
  }
  data = do.call("rbind",data)
  
  results = list(fileIds=matrix(1:(fileCounter-1),ncol=length(conditions),dimnames=list(c(),conditions)),fileNames=fileNames,fileChannelNames=fileChannelNames,fileReagentNames=fileReagentNames)
  
  if (!is.null(transformColumns)){
    if (any(!is.numeric(transformColumns))){
      containedCols = setdiff(transformColumns,colnames(data))
      if (length(containedCols)>0){
        stop(paste("Transform cols",paste(containedCols,collapse=", "),"not found. Valid channel names:",paste(colnames(data),collapse=", ")))
      }  
    }
    data[,transformColumns] = asinh(data[,transformColumns]/transformCofactor);
    results$transformColumns = transformColumns
    results$transformCofactor = transformCofactor
  }
  
  if (!is.null(scaleColumns)){
    if (any(!is.numeric(scaleColumns))){
      containedCols = setdiff(scaleColumns,colnames(data))
      if (length(containedCols)>0){
        stop(paste("Scale cols",paste(containedCols,collapse=", "),"not found. Valid channel names:",paste(colnames(data),collapse=", ")))
      }  
    }
    results$scaleColumns = scaleColumns
    results$scaleColumns.mean = apply(data[,scaleColumns],2,mean)  
    results$scaleColumns.SD = apply(data[,scaleColumns],2,sd)  
    data[,scaleColumns] = apply(data[,scaleColumns],2,scale)  
  }
  
  results$data = data
  

  class(results) = "citrus.combinedFCSSet"
  
  return(results);
}


#' @export
#' @name citrus.readFCSSet
print.citrus.combinedFCSSet = function(citrus.combinedFCSSet){
  cat(paste0("Number Of Files: ",length(citrus.combinedFCSSet$fileNames),"\n"))
  cat(paste0("Number Of Events: ",nrow(citrus.combinedFCSSet$data),"\n"))
  cat(paste0("Parameter Names:\n",paste0("\t",colnames(citrus.combinedFCSSet$data),collapse="\n")))
}

#' Masks a citrus.combinedFCSSet
#' 
#' Masks a citrus.combinedFCSset to include data from a subset of file ids. 
#' 
#' @param citrus.combinedFCSSet A citrus.combinedFCSSet object.
#' @param fileIds Vector of file IDs for which to retain data. 
#' 
#' @return A \code{citrus.combinedFCSSet} object.
#' @author Robert Bruggner
#' @export
#' 
#' @seealso \code{\link{citrus.readFCSSet}}
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read citrus.combinedFCSSet
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # Mask
#' maskedFCSSet = citrus.maskCombinedFCSSet(citrus.combinedFCSSet,fileIds=1:3)
#' 
#' # Check - should be 1,2,3
#' unique(maskedFCSSet$data[,"fileId"])
citrus.maskCombinedFCSSet = function(citrus.combinedFCSSet,fileIds){

  # Keep any row that has an included file
  keepFileIds = citrus.combinedFCSSet$fileIds[unique(which(apply(citrus.combinedFCSSet$fileIds,2,"%in%",fileIds),arr.ind=T)[,1]),,drop=F]
  citrus.combinedFCSSet$data = do.call("rbind",lapply(fileIds,subsetByFile,data=citrus.combinedFCSSet$data))
  citrus.combinedFCSSet$fileIds = citrus.combinedFCSSet$fileIds[unique(which(apply(citrus.combinedFCSSet$fileIds,2,"%in%",fileIds),arr.ind=T)[,1]),,drop=F]
  
  return(citrus.combinedFCSSet)
}

subsetByFile = function(data,fileId){
  subset(data,data[,"fileId"]==fileId)
}

#' Select clusters for further analysis
#' 
#' Selects clusters from a set of clusters identified by a clustering method investigate for stratifying signal.
#' 
#' @param citrus.clustering A \code{citrus.clustering} object.
#' @param method Method for determining which clusters from a clustering should be analyzed for stratifying signal.
#' @param ... Other parameters passed to specific cluster selection methods.
#' 
#' @return A vector of cluster IDs. 
#' 
#' @author Robert Bruggner
#' @export
#' @seealso \code{\link{citrus.selectClusters.minimumClusterSize}}.
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
#' # Select clusters that contain at least 1% of clustered events.
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering,minimumClusterSizePercent=0.01)
citrus.selectClusters = function(citrus.clustering,method="minimumClusterSize",...){
  do.call(paste0("citrus.selectClusters.",method),args=list(citrus.clustering=citrus.clustering,...))
}

#' Selects clusters for endpoint analysis
#' 
#' Selects clusters for endpoint analysis by cluster size. Selected clusters must have a minimum number of 
#' cells in them as a proportion of the total number of clustered events. If \eqn{n} total events are clustered,
#' clusters contatining at least \eqn{n * minimumClusterSizePercent} events are selected.
#' 
#' @param citrus.clustering A \code{citrus.clustering} object.
#' @param minimumClusterSizePercent The percentage \eqn{0 < x < 1} of the total number of clustered events a cluster 
#' must contain in order to be selected. 
#' @param ... Other arguments (ignored).
#' 
#' @author Robert Bruggner
#' @export
#' @seealso \code{\link{citrus.selectClusters}}
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
#' # Select clusters that contain at least 1% of clustered events.
#' largeEnoughClusters = citrus.selectClusters.minimumClusterSize(citrus.clustering,minimumClusterSizePercent=0.01)
citrus.selectClusters.minimumClusterSize =function(citrus.clustering,minimumClusterSizePercent=0.05,...){
  clusterSizes = sapply(citrus.clustering$clusterMembership,length)
  minimumClusterSize = (length(citrus.clustering$clusterMembership)+1)*minimumClusterSizePercent
  return(which(clusterSizes>=minimumClusterSize))
}


citrus.convertConditionMatrix = function(conditionMatrix){
  conditions = list();
  for (i in 1:nrow(conditionMatrix)){
    for (j in 1:ncol(conditionMatrix)){
      if (conditionMatrix[i,j]){
        conditions = append(conditions,list(unique(c(rownames(conditionMatrix)[i],colnames(conditionMatrix)[j]))))
      }
    }
  }
  return(conditions)
}

#' Get Citrus package version
#' 
#' Get Citrus package version.
#' 
#' @return Citrus package version
#' @author Robert Bruggner
#' @export
citrus.version = function(){
  as.character(packageVersion("citrus"))
}

#' Counts FCS File events
#' 
#' Counts FCS events in all FCS files in a data directory
#' 
#' @param dataDirectory Directory from which to read FCS files. 
#' @param ... Other parameters passed to citrus.readFCS
#' 
#' @author Robet Bruggner
#' @export
#' 
#' @return A matrix detailing the number of events and parameters in each directory FCS file.
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Get FCS Event Counts
#' citrus.fileEventCount(dataDirectory)
citrus.fileEventCount = function(dataDirectory,...){
  lengths = list();
  
  for (fcsFile in list.files(dataDirectory,pattern=".fcs",ignore.case=T)){
    print(paste("Reading",fcsFile))
    lengths[[fcsFile]] = suppressWarnings(dim(citrus.readFCS(file.path(dataDirectory,fcsFile),...)))
  }
  return(do.call("rbind",lengths))
}

#' List of computable feature types
#' 
#' Returns valid types of descriptive features that Citrus can compute
#' 
#' @return Vector of feature types Citrus can compute.
#' 
#' @author Robert Bruggner
#' @export
citrus.featureTypes = function(){
  return(c("abundances","medians"))
}

#' List possible model families 
#' 
#' Returns valid model family types that Citrus can compute.
#' 
#' @return Vector of model families that Citrus can compute.
#' 
#' @author Robert Bruggner
#' @export
citrus.familyList = function(){
  return(c("classification","quantiative"))
}


#' List possible model types
#' 
#' Returns valid model types that Citrus can compute.
#' 
#' @return Vector of model types that Citrus can compute.
#' 
#' @author Robert Bruggner
#' @export
citrus.modelTypes = function(){
  return(c("pamr","glmnet","sam"))
}

#' Read an FCS file 
#' 
#' Reads an FCS file while suppressing warnings.
#' 
#' @param filePath The full path to the FCS file to be read.
#' @param ... Other arguments to be passed to read.FCS.

#' @return a \code{\link{flowCore}} flowFrame object.
#' @author Robert Bruggner
#' @export 
#' @seealso \code{\link{flowCore}}
citrus.readFCS = function(filePath,...){
  addtlArgs = list(...)
  dataset=1
  if ("dataset" %in% names(addtlArgs))
    dataset=addtlArgs[["dataset"]]
  
  which.lines=NULL
  if ("which.lines" %in% names(addtlArgs))
    which.lines=addtlArgs[["which.lines"]]
  
  fcs = tryCatch({
    suppressWarnings(read.FCS(filePath,dataset=dataset,which.lines=which.lines))
  }, error = function(e){
    if (grepl("Please set argument 'emptyValue' as",e$message)){
      suppressWarnings(read.FCS(filePath,dataset=dataset,which.lines=which.lines,emptyValue=F))
    } else {
      stop(e$message)
    }
  }) 
  return(fcs)
}


#' Exports cluster events to file
#' 
#' Exports cluster events from each file
#' 
#' @param clusterId ID of cluster to be exported
#' @param citrus.clustering \code{citrus.clustering} object to export clusters from.
#' @param citrus.combinedFCSSet \code{citrus.combinedFCSSet} object that contains data to be exported.
#' @param outputDirectory Where cluster data should be exported to.
#' @param conditions If not \code{NULL}, only export data from the specified conditions.
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
#' # Export Cluster
#' # citrus.exportCluster(clusterId=19500,citrus.clustering,citrus.combinedFCSSset,outputDirectory="/PATH/TO/DIRECTORY/")
citrus.exportCluster = function(clusterId,citrus.clustering,citrus.combinedFCSSet,outputDirectory,conditions=NULL){
  
  # Get condition file ids for export
  if (is.null(conditions)){
    conditions = colnames(citrus.combinedFCSSet$fileIds)
  }
  fileIds = as.vector(citrus.combinedFCSSet$fileIds[,conditions])
  
  # Get data from all files within cluster
  clusterData = citrus.combinedFCSSet$data[citrus.clustering$clusterMembership[[clusterId]],]
  
  # Data should be unscaled before export, if scaled at read time.
  if (!is.null(citrus.combinedFCSSet$scaleColumns)){
    clusterData[,citrus.combinedFCSSet$scaleColumns] = t((t(clusterData[,citrus.combinedFCSSet$scaleColumns])*citrus.combinedFCSSet$scaleColumns.SD)+citrus.combinedFCSSet$scaleColumns.mean)
  }
  
  # Data should be untransformed before export, if it was transformed at read time.
  if (!is.null(citrus.combinedFCSSet$transformColumns)){
    clusterData[,citrus.combinedFCSSet$transformColumns] = sinh(clusterData[,citrus.combinedFCSSet$transformColumns])*citrus.combinedFCSSet$transformCofactor
  }
  
  # Export cluster data from each file at a time
  for (fileId in fileIds){
    # Get data from cluster in file
    clusterFileData = clusterData[clusterData[,"fileId"]==fileId,,drop=F]
    fileName = paste0(citrus.combinedFCSSet$fileNames[fileId],".cluster_",clusterId,".fcs")
    
    # TODO: Fix to include appropriate chanel names and descriptions instead of just one or the other
    write.FCS(flowFrame(clusterFileData),filename=file.path(outputDirectory,fileName))
  }
  
}

#' UI for exporting file cluster events 
#' 
#' UI for exporting file cluster events. When run, a selection dialog will prompt the user to select
#' a file. The user should select a "citrusClustering.rData" file from a citrusOutput directory.
#' The user should then enter a comma-separated list of clusters IDs to be exported. A directory
#' called 'clusterExportData' will be created in the same directory as the selected clustering file
#' and exported clusters will be placed there. 
#' 
#' @author Robert Bruggner
#' @export
citrus.exportClusterUI = function(){
  # Select the clustering to export
  clusteringOutputPath = file.choose()
  cat(paste0("Opening Clustering: ",clusteringOutputPath))
  load(clusteringOutputPath);
  
  # Get cluster ids to export from 
  clusterExportIdsString = readline("Enter comma-separated list of cluster ids for export:")
  exportClusterIds = as.numeric(unlist(strsplit(clusterExportIdsString,",")))
  
  # Create directory to put output
  exportDir = file.path(dirname(clusteringOutputPath),"clusterExportData")
  dir.create(exportDir,showWarnings = F)
  
  # Export clusters
  sapply(exportClusterIds,citrus.exportCluster,citrus.clustering=citrus.foldClustering$allClustering,citrus.combinedFCSSet=citrus.combinedFCSSet,outputDirectory=exportDir)
}

#' Graphical Interface to checking file parameter consistency.
#' 
#' Graphical Interface to checking file parameter consistency.
#' 
#' @details When running \code{citrus.checkFileParameterConsistencyUI()}, the user will first be prompted to select a file.
#' The user should then select a single FCS file from the directory where all FCS files reside. Citrus will then check the 
#' consistency of parameter names and counts for all FCS files in that directory.
#' 
#' @author Robert Bruggner
#' @export
citrus.checkFileParameterConsistencyUI = function(){
  fcsFilePath = file.choose()  
  dataDirectory = dirname(fcsFilePath)
  citrus.checkFileParameterConsistency(dataDirectory)
}

#' Checks consistency of FCS File parameters.
#' 
#' Checks consistency of FCS File parameter names. Will report if FCS files have differing numbers of
#' parameters or parameters have different names.
#' 
#' @param dataDirectory Full path to directory where FCS files are. 
#' 
#' @return Returns \code{0} if file parameters are consistent, \code{-1} otherwise.
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Check parameter consistency. 
#'  citrus.checkFileParameterConsistency(dataDirectory)
citrus.checkFileParameterConsistency = function(dataDirectory){
  
  #Get list of files in data directory
  fcsFiles = list.files(dataDirectory,pattern=".FCS",ignore.case = T)
  fileParameters = lapply(fcsFiles,citrus.getFileParameters,dataDirectory)
  names(fileParameters)=fcsFiles
  fileParameters 
  
  # Check for unequal numbers of parameters in files 
  if (length(unique(sapply(fileParameters,length)))>1){
    print("One or more FCS files has a different number of channels. Please fix before proceeding.")
    print("Number of parameters per file:")
    print(sapply(fileParameters,length))
    return(-1)
  }
  
  # Check for unlike names in fcs parameters names
  if (any(sapply(apply(do.call("rbind",fileParameters),2,unique),length)>1)){
    print("One or more FCS files has different channel names. Please fix before proceeding.")
    
    problemChannels = which(sapply(apply(do.call("rbind",fileParameters),2,unique),length)>1)
    for (problemChannel in names(problemChannels)){
      print(paste0("Channel Name: ",problemChannel))
      print( 
        paste0("Unique Names: ",paste(apply(do.call("rbind",fileParameters),2,unique)[[problemChannel]],collapse=", "))
      )
    }
    return(-1)
  }
  
  print("FCS File parameters consistent.")
  return(0)
}
