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
#' @param ... Other undocumented parameters.
#' 
#' @return An object of type \code{citrus.combinedFCSSet}.
#' \item{data}{Combined FCS measurements from files listed in fileList, with paraemters \code{fileEventNumber} and \code{fileId} added. The former
#' reports the event number in the original source file and the latter records the file from which the event was sampled.}
#' \item{fileIds}{A matrix recording the assigned fileIds for each file in the fileList argument.}
#' \item{fileNames}{The corresponding file names for each file, indexed by fileId.}
#' \item{fileChannelNames}{Names of parameters in each read file.}
#' \item{fileReagentNames}{Description of parameters in each read file.}
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
citrus.readFCSSet = function(dataDirectory,fileList,fileSampleSize=NULL,transformColumns=NULL,transformCofactor=5,scaleColumns=NULL,useChannelDescriptions=F,...){
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
      
      fileChannelNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$name)
      fileReagentNames[[conditions[i]]][[fileName]]=as.vector(pData(parameters(fcsFile))$desc)
      fcsData = cbind(fcsData,fileEventNumber=1:nrow(fcsData),fileId=fileCounter);
      fileCounter=fileCounter+1;
      
      if (!is.null(transformColumns)){
        if (any(!is.numeric(transformColumns))){
          containedCols = setdiff(transformColumns,colnames(fcsData))
          if (length(containedCols)>0){
            stop(paste("Transform cols",paste(containedCols,collapse=", "),"not found. Valid channel names:",paste(colnames(fcsData),collapse=", ")))
          }  
        }
        fcsData[,transformColumns] = asinh(fcsData[,transformColumns]/transformCofactor);
      }
      
      if ((!is.null(fileSampleSize))&&(fileSampleSize<nrow(fcsData))){
        cat(paste("\tSampling",fileSampleSize,"events.\n"))
        fcsData = fcsData[sort(sample(1:nrow(fcsData),fileSampleSize)),] 
      }
      
      if ("scaleSampleCols" %in% names(addtlArgs)){
        fcsData[,addtlArgs[["scaleSampleCols"]]] = apply(fcsData[,addtlArgs[["scaleSampleCols"]]],2,scale,center=F)  
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
  
  if (!is.null(scaleColumns)){
    fcsData[,scaleColumns] = apply(fcsData[,scaleColumns],2,scale)  
  }
  
  results = list(data=data,fileIds=matrix(1:(fileCounter-1),ncol=length(conditions),dimnames=list(c(),conditions)),fileNames=fileNames,fileChannelNames=fileChannelNames,fileReagentNames=fileReagentNames)
  class(results) = "citrus.combinedFCSSet"
  
  return(results);
}

#' Prints a summary of a citrus.combinedFCSSet.
#' 
#' @param citrus.combinedFCSSet A citrus.combinedFCSSet object.
#' 
#' @method print citrus.combinedFCSSet
#' @S3method print citrus.combinedFCSSet
#' 
#' @author Robert Bruggner
#' @export
print.citrus.combinedFCSSet = function(citrus.combinedFCSSet,...){
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
  
  results = list(data=do.call("rbind",lapply(fileIds,.subsetByFile,data=citrus.combinedFCSSet$data)),
             fileIds=citrus.combinedFCSSet$fileIds[unique(which(apply(citrus.combinedFCSSet$fileIds,2,"%in%",fileIds),arr.ind=T)[,1]),,drop=F],
             fileNames=citrus.combinedFCSSet$fileNames,
             fileChannelNames=citrus.combinedFCSSet$fileChannelNames,
             fileReagentNames=citrus.combinedFCSSet$fileReagentNames,
             call = match.call() 
        )
  class(results) = "citrus.combinedFCSSet"
  return(results)
}

.subsetByFile = function(data,fileId){
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
#' cells in them as a proportion of the total number of clustered events. If \code{n} total events are clustered,
#' clusters contatining at least \code{n * minimumClusterSizePercent} events are selected.
#' 
#' @param citrus.clustering A \code{citrus.clustering} object.
#' @param minimumClusterSizePercent The percentage (0 < x < 1) of the total number of clustered events a cluster 
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


#citrus.assembleHandGates = function(dataDir,filePopulationList,conditionComparaMatrix=NULL,fileSampleSize=NULL,transformColumns=NULL,transformCofactor=5){
  
#  if (!is.null(conditionComparaMatrix)){
#    allConditions = citrus.convertConditionMatrix(conditionComparaMatrix) 
#    # Perform Internal Consistency Checks!
#  } else {
#    allConditions = names(filePopulationList)
#  }
  
#  results = list();

#  for (conditions in allConditions){
    
#    cda = data.frame();
#    eventCounter=0;
#    populationCounter=1;  
#    fileNames=list()  
#    clusterAssignments = list();
#   fileId=1;  
#    fileChannelNames = list()
#    fileReagentNames = list()
#    fileIds = list();
#    populationClusterMap = list();
#    for (condition in conditions){
#      
#      nPatients = nrow(filePopulationList[[condition]])
#      baseFileId = (match(condition,conditions)-1)*nPatients
#      fileIds[[condition]] = (baseFileId+1):(baseFileId+nPatients)
#      for (populationName in colnames(filePopulationList[[condition]])){
#        populationEvents = c();
#        for (patientId in 1:nPatients){
#          fcsFileName = filePopulationList[[condition]][patientId,populationName]
#          fileId = baseFileId+patientId
#          print(paste("Reading",fcsFileName))
#          fcsFile = citrus.readFCS(file.path(dataDir,fcsFileName),...)
#          fileChannelNames[[fcsFileName]] = colnames(fcsFile)
#          fileReagentNames[[fcsFileName]] = pData(parameters(fcsFile))$desc
#          fcsFileData = exprs(fcsFile)
#          if (!is.null(fileSampleSize)){
#            if (fileSampleSize < nrow(fcsFileData)){
#              fcsFileData = fcsFileData[sample(1:nrow(fcsFileData),fileSampleSize),]
#            }
#          }
#          lastEvent = nrow(cda)
#          absoluteEventIds = (lastEvent+1):(eventCounter+nrow(fcsFileData))
#          cda = rbind(cda,data.frame(fcsFileData,fileId=fileId))
#          populationEvents = c(populationEvents,absoluteEventIds)
#          eventCounter = eventCounter+nrow(fcsFileData)
#          fileNames[[fileId]]=fcsFileName
#          fileId=fileId+1;
#        }
#        
#        
#        if (populationName %in% names(populationClusterMap)){
#          clusterAssignments[[populationClusterMap[[populationName]]]] = c(clusterAssignments[[populationClusterMap[[populationName]]]],populationEvents)
#        } else {
#          clusterAssignments[[populationCounter]]=populationEvents;
#          populationClusterMap[[populationName]]=populationCounter
#          populationCounter = populationCounter+1;
#        }
#        
#      }
#    }
#    
#    if (!is.null(transformColumns)){
#      cda[,transformColumns] = cda[,transformColumns]/transformCofactor
#    }
#    
#    results[[paste(conditions,collapse="_vs_")]] = list(folds="all",
#                                                        foldsClusterAssignments=list(all=clusterAssignments),
#                                                        conditions=conditions,
#                                                        citrus.dataArray=list(data=as.matrix(cda),
#                                                                              fileNames=unlist(fileNames,use.names=F),
#                                                                              fileIds=do.call("cbind",fileIds),
#                                                                              fileChannelNames=fileChannelNames,
#                                                                              fileReagentNames=fileReagentNames),
#                                                        clusterPopulationNames=names(populationClusterMap));
#  }
#  return(results)
#}

#citrus.leavoutFold = function(x,y,leaveoutSize){
#  groupCounts = table(y)/length(y)*leaveoutSize
#  leaveout = list()
#  for (group in unique(y)){
#    leaveout[[group]] = sample(which(y==group),groupCounts[which(names(groupCounts)==group)])
#  }
#  return(as.vector(unlist(leaveout)))
#}

citrus.formatDecimal = function(x){
  sprintf("%1.2f", x)
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

citrus.version = function(){
  as.character(packageVersion("citrus"))
}

citrus.fileEventCount = function(dataDir,...){
  lengths = list();
  
  for (fcsFile in list.files(dataDir,pattern=".fcs",ignore.case=T)){
    print(paste("Reading",fcsFile))
    lengths[[fcsFile]] = suppressWarnings(dim(citrus.readFCS(file.path(dataDir,fcsFile),...)))
  }
  return(do.call("rbind",lengths))
}

citrus.featureTypes = function(){
  return(c("abundances","medians"))
}

citrus.familyList = function(){
  return(c("classification","survival","quantiative"))
}

citrus.modelTypes = function(){
  return(c("pamr","glmnet","sam"))
}

#' Read an FCS file 
#' 
#' \code{citrus.readFCS} reads an FCS file while suppressing warnings
#' 
#' @param filePath the full path to the FCS file to be read.
#' @param ... other arguments to be passed to read.FCS
#' @return a flowCore flowFrame object
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

citrus.exportConditionClusters = function(conditionClusterIds,preclusterResult,outputDir,sampleIds=NULL){
  for (conditionName in names(conditionClusterIds)){
    for (clusterId in conditionClusterIds[[conditionName]]){
      nFolds = length(preclusterResult[[conditionName]]$foldsCluster)
      outputFile = file.path(outputDir,paste0(conditionName,"-cluster_",clusterId,".fcs"))
      citrus.exportCluster(clusterId,
                           data=preclusterResult[[conditionName]]$citrus.dataArray$data,
                           clusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments[[nFolds]],
                           outputFile=outputFile,
                           sampleIds=sampleIds)
    }
  }
}

#' @title Exports a cluster of cells to an FCS file 
#' 
#' @description Export cells in a cluster to an FCS file
#' 
#' @param clusterId The numeric ID of the cluster to be exported
#' @param data The matrix of combined data that was clustered
#' @param clusterAssignments A list with each element containin the indices of cells in each cluster. Generated by citrus.
citrus.exportCluster = function(clusterId,data,clusterAssignments,outputFile,sampleIds=NULL){
  clusterData = data[clusterAssignments[[clusterId]],]
  if (!is.null(sampleIds)){
    clusterData = clusterData[clusterData[,"fileId"]%in%sampleIds,]
  }
  if (nrow(clusterData)==0){
    warning(paste("No data for cluster ",clusterId,"in selected files. Not writing to disk."))
  } else {
    write.FCS(x=flowFrame(exprs=clusterData),filename=outputFile)  
  } 
}
