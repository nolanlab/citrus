.extractConditionLargeEnoughClusters = function(condition,foldFeatures){
  foldFeatures[[condition]]$foldLargeEnoughClusters[[which(foldFeatures[[condition]]$folds=="all")]]
}

#citrus.buildFeatures = function(preclusterResult,outputDir,featureTypes=c("abundances"),minimumClusterSizePercent=0.05,largeEnoughClusters=NULL,conditionCluster=NULL,...){  
  
#  addtlArgs = list(...)
  
  # Error check before we actually start the work.
#  if ((!all(featureTypes %in% citrus.featureTypes()))||(length(featureTypes)<1)){
#    stop(paste("featureTypes must be 1 or more of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
#  }
#  
#  if (("medians" %in% featureTypes)&&(!("medianColumns" %in% names(list(...))))){
#    stop("medianColumns argument must be specified to calculate cluster medians.")
#  }
#  
#  if (("emDists" %in% featureTypes)&&(!("emdColumns" %in% names(list(...))))){
#    stop("emdColumns argument must be specified to calculate cluster emDists.")
#  }
#  
#  
#  if (!file.exists(outputDir)){
#    stop(paste("Output directory",outputDir,"not found."))
#  }
#  
#  if (!is.null(conditionCluster)){
#    featureRes = parLapply(conditionCluster,names(preclusterResult),citrus.buildConditionFeatures,preclusterResult=preclusterResult,featureTypes=featureTypes,minimumClusterSizePercent=minimumClusterSizePercent,largeEnoughClusters=largeEnoughClusters,outputDir=outputDir,...)      
#  } else {
#    featureRes = lapply(names(preclusterResult),citrus.buildConditionFeatures,preclusterResult=preclusterResult,featureTypes=featureTypes,minimumClusterSizePercent=minimumClusterSizePercent,largeEnoughClusters=largeEnoughClusters,outputDir=outputDir,...)  
#  }
#  
# names(featureRes) = names(preclusterResult)
#  return(featureRes)
#}

citrus.buildFeatures = function(citrus.combinedFCSSet,citrus.clustering,featureType="abundances",minimumClusterSizePercent=0.05,...){  

  addtlArgs = list(...)

  #Error check before we actually start the work.
  if ((!all(featureType %in% citrus.featureTypes()))||(length(featureType)<1)){
    stop(paste("featureType must be one of the following:",paste(citrus.featureTypes(),collapse=", "),"\n"))
  }
  
  if (("medians" %in% featureType)&&(!("medianColumns" %in% names(list(...))))){
    stop("medianColumns argument must be specified to calculate cluster medians.")
  }
  
  clusterSizes = sapply(citrus.clustering$clusterMembership,length)
  largeEnoughClusters = which(clusterSizes>=nrow(citrus.combinedFCSSet$data)*minimumClusterSizePercent)
  do.call(paste0("citrus.calculateFeature.",featureType),args=list(clusterIds=largeEnoughClusters,citrus.clustering=citrus.clustering,citrus.combinedFCSSet=citrus.combinedFCSSet,...))
}

################################
# Abundance Features
################################
citrus.calculateFeature.abundances = function(clusterIds,citrus.clustering,citrus.combinedFCSSet,...){
  eventFileIds = citrus.combinedFCSSet$data[,"fileId"]
  fileIds = unique(eventFileIds)
  citrus.calculateFileClusterAbundance(clusterId=clusterIds[1],fileId=fileIds[1],clusterAssignments=citrus.clustering$clusterMembership,eventFileIds=eventFileIds)
  res = mcmapply(citrus.calculateFileClusterAbundance,
          clusterId=rep(clusterIds,length(fileIds)),
          fileId=rep(fileIds,each=length(clusterIds)),
          MoreArgs=list(
                        clusterAssignments=citrus.clustering$clusterMembership,
                        eventFileIds=eventFileIds)
          ,...)
  return(
    matrix(res,ncol=length(clusterIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames,paste("cluster",clusterIds,"abundance")))
  )
}

citrus.calculateFileClusterAbundance = function(clusterId,fileId,clusterAssignments,eventFileIds){
  sum(which(eventFileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((eventFileIds==fileId))
}

################################
# Median Features
################################
citrus.calculate.calculateFeature.medians = function(clusterIds,citrus.clustering,citrus.combinedFCSSet,...){
  if !("medianColumns" %in% names(list(...)))
    stop("medianColumns argument missing")
  
  medianColumns = list(...)[["medianColumns"]]
  fileIds = unique(citrus.combinedFCSSet$data[,"fileId"])
  res = mcmapply(citrus.calculateFileClusterMedian,
                 clusterId=rep(rep(clusterIds,length(medianColumns)),length(fileIds)),
                 fileId=rep(fileIds,each=(length(clusterIds)*length(medianColumns))),
                 medianColumn=rep(rep(medianColumns,each=length(clusterIds)),length(fileIds)),
                 MoreArgs=list(
                   data=citrus.combinedFCSSet$data,
                   clusterAssignments=citrus.clustering$clusterMembership)
        ,...)
  
  return(
    matrix(res,nrow=length(fileIds),byrow=T,dimnames=list(citrus.combinedFCSSet$fileNames,paste("cluster",rep(clusterIds,length(medianColumns)),rep(medianColumns,each=length(clusterIds)),"median")))
  )
}

citrus.calculateFileClusterMedian = function(clusterId,fileId,medianColumn,data,clusterAssignments){
  clusterData = data[clusterAssignments[[clusterId]],]
  clusterFileDataValues = clusterData[clusterData[,"fileId"]==fileId,medianColumn]
  if (length(clusterFileDataValues)<3){
    return(0)
  } else {
    return(median(clusterFileDataValues))
  }
}
  

#citrus.calculateFeature.medians = function(foldsFileIds,clusterIds,clusterAssignments,data,citrus.dataArray,...){
#  addtlArgs = list(...)
#  if ("mc.cores" %in% names(addtlArgs)){
#    features = do.call("rbind",mclapply(foldsFileIds,citrus.calculateFileClustersMedians,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))    
#  } else {
#    features = t(sapply(foldsFileIds,citrus.calculateFileClustersMedians,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))    
#  }
#  
#  rownames(features) = citrus.dataArray$fileNames[foldsFileIds]
# return(features)
#}

#citrus.calculateFileClustersMedians = function(fileId,clusterIds,clusterAssignments,data,...){
#  addtlArgs = list(...)
#  if (!("medianColumns" %in% names(addtlArgs))){
#    stop("medianColumns argument must be specified to compute cluster medians.")
#  }
  
#  if (any(!is.numeric(addtlArgs[["medianColumns"]]))){
#    containedCols = setdiff(addtlArgs[["medianColumns"]],colnames(data))
#    if (length(containedCols)>0){
#      stop(paste("Cluster cols",paste(containedCols,collapse=", "),"not found. Valid channel names:",paste(colnames(data),collapse=", ")))
#    }  
#  }
#  unlist(lapply(clusterIds,citrus.calculateFileClusterMedian,clusterAssignments=clusterAssignments,fileId=fileId,data=data,medianColumns=addtlArgs[["medianColumns"]]))
#}

#citrus.calculateFileClusterMedian = function(clusterId,clusterAssignments,fileId,data,medianColumns){
#  include = data[clusterAssignments[[clusterId]],]
  
  # Are there zero cells assigned to our cluster?
#  if (length(clusterAssignments[[clusterId]])==0){
#    medians = rep(0,length(medianColumns))
  
#  } else if (length(clusterAssignments[[clusterId]])==1) {
    # Are there 1 cells assigned to our cluster?
    
#    if (fileId==include["fileId"]){
#      medians = rep(0,length(medianColumns))
#    } else {
#      medians = include[medianColumns]
#    }
#    
#  } else {
#    # Are there more than 1 cells assigned to our cluster?
#    
#    include = include[include[,"fileId"]==fileId,medianColumns]
#    if (length(medianColumns)>1){
#      # do we have data?
#      if (length(include)>0){
#        if (is.null(nrow(include))){
#          medians = include  
#        } else {
#          medians = apply(include,2,median) 
#        }
#      } else {
#        medians = (rep(0,length(medianColumns)))
#      }
#    } else if (length(medianColumns)==1){
#      if (length(include)>0){
#        medians = median(include)
#     } else {
#        medians=0
#      }
#    }
#  }
#  
#  names(medians) = paste(paste(paste("cluster",clusterId),medianColumns),"median")
#  return(medians)
#}


citrus.calculateConditionFeatureDifferences = function(data,condition1,condition2){
  x = data[condition2,]-data[condition1,]
  rownames(x)=apply(cbind(rownames(data[condition2,]),rownames(data[condition1,])),1,paste,collapse=" - ")
  return(x)
}

citrus.calculateFoldLargeEnoughClusters = function(index,foldsClusterAssignments,folds,citrus.dataArray,conditions,minimumClusterSizePercent=0.05){
  if ((length(folds[[index]])==1) && (folds[[index]]=="all")){
    foldsDataLength = sum(citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions])
  } else {
    foldsDataLength = sum(citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[-folds[[index]],conditions])
  }
  minimumClusterSize = floor(minimumClusterSizePercent*foldsDataLength)
  return(citrus.calculateLargeEnoughClusters(foldsClusterAssignments[[index]],minimumClusterSize))
}
  
citrus.calculateLargeEnoughClusters = function(clusters,minimumClusterSize){
  clusterLengths = lapply(clusters,length)
  return(which(clusterLengths >= minimumClusterSize))
}

citrus.getNonOverlappingFeatures = function(index,foldFeatures,leftoutFeatures){
  c(setdiff(colnames(foldFeatures[[index]]),colnames(leftoutFeatures[[index]])),setdiff(colnames(leftoutFeatures[[index]]),colnames(foldFeatures[[index]])))
}

citrus.removeFeatures = function(index,foldFeatures,nonOverlappingFeatures){
  ff = foldFeatures[[index]]
  remove = which(colnames(ff)%in%nonOverlappingFeatures[[index]])
  if (length(remove)>0){
    ff = ff[,-remove]
  }
  return(ff)
}


citrus.getCVMinima = function(modelType,thresholdCVRates,fdrRate=0.01){
  cvPoints=list();
  if (modelType=="sam"){
    cvPoints[["fdr_0.10"]]=10
    cvPoints[["fdr_0.05"]]=5
    cvPoints[["fdr_0.01"]]=1
  } else {
    errorRates = thresholdCVRates[[modelType]]$cvm
    SEMs = thresholdCVRates[[modelType]]$cvsd
    FDRRates = thresholdCVRates[[modelType]]$fdr
    cvPoints[["cv.min.index"]] = min(which(errorRates==min(errorRates,na.rm=T)))
    cvPoints[["cv.min"]] = thresholdCVRates[[modelType]]$threshold[cvPoints[["cv.min.index"]]]
    cvPoints[["cv.1se.index"]] = min(which(errorRates<=(errorRates[cvPoints[["cv.min.index"]]]+SEMs[cvPoints[["cv.min.index"]]])))
    cvPoints[["cv.1se"]] = thresholdCVRates[[modelType]]$threshold[cvPoints[["cv.1se.index"]]]
    if (!is.null(FDRRates)) {
      if (any(FDRRates<fdrRate)){
        if (length(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates,na.rm=T))))>0){
          cvPoints[["cv.fdr.constrained.index"]] = max(intersect(which(FDRRates<0.01),which(errorRates==min(errorRates,na.rm=T))))
          cvPoints[["cv.fdr.constrained"]] = thresholdCVRates[[modelType]]$threshold[cvPoints[["cv.fdr.constrained.index"]]]
        }
      }
      
    }
  }
  return(cvPoints)
}


citrus.extractModelFeatures = function(modelType,cvMinima,foldModels,foldFeatures,regularizationThresholds,family){
  res = list();
  nAllFolds = length(foldModels[[modelType]])
  finalModel = foldModels[[modelType]][[nAllFolds]]
  for (cvPoint in names(cvMinima[[modelType]])[!grepl("index",names(cvMinima[[modelType]]))]){
    threshold = cvMinima[[modelType]][[cvPoint]]
    thresholdIndex = cvMinima[[modelType]][[paste(cvPoint,"index",sep=".")]]
    if (modelType=="pamr"){
      if (finalModel$nonzero[thresholdIndex]>0){
        f = pamr.listgenes(fit=finalModel,data=list(x=t(foldFeatures[[nAllFolds]]),geneids=colnames(foldFeatures[[nAllFolds]])),threshold=threshold)  
        f = as.vector(f[,1])
        res[[cvPoint]][["features"]] = f
        res[[cvPoint]][["clusters"]] = sort(unique(as.numeric(do.call("rbind",strsplit(f,split=" "))[,2])))  
      } else {
        res[[cvPoint]][["features"]] = NULL
        res[[cvPoint]][["clusters"]] = NULL
      }
      
    } else if (modelType=="glmnet"){
      # THIS NEEDS TO BE FIXED IN ORDER TO SUPPORT MULTINOMIAL REGRESSION WITH GLMNET
      f = as.matrix(predict(finalModel,newx=foldFeatures[[nAllFolds]],type="coefficient",s=threshold))
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


#citrus.buildConditionFeatures = function(conditionName,preclusterResult,featureTypes,minimumClusterSizePercent,largeEnoughClusters,outputDir,...){
#  cat(paste("Building features for condition",conditionName,"\n"))
#  conditions=preclusterResult[[conditionName]]$conditions

#  folds = preclusterResult[[conditionName]]$folds
#  nAllFolds=length(folds)
#  nFolds=nAllFolds-1

#  if (is.null(largeEnoughClusters)){
#    cat("Calculating Fold Large Enough Clusters\n")
#    foldLargeEnoughClusters = lapply(1:nAllFolds,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,minimumClusterSizePercent=minimumClusterSizePercent)  
#  } else {
#    foldLargeEnoughClusters = list(mappingResult=largeEnoughClusters[[conditionName]])
#  }

#  cat("Calculating Features\n")
#  foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,...)
#  if (any(do.call("c",lapply(lapply(foldFeatures,dim),is.null)))){
#    stop("No Features Calculated.")
#  }
#foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,emdColumns=emdColumns)

#Normalize features... Sometimes EMD's aren't calculated. Need a better way to handle this.
#  if (!(folds[[1]][1] %in% c("all","mappingResult"))){
#    leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=preclusterResult[[conditionName]]$citrus.dataArray,foldsClusterAssignments=preclusterResult[[conditionName]]$leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,...)
#leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=featureTypes,folds=folds,citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,calculateLeaveoutData=T,emdColumns=emdColumns)

#    nof = lapply(1:nFolds,citrus.getNonOverlappingFeatures,foldFeatures=foldFeatures,leftoutFeatures=leftoutFeatures)
#    foldFeatures[1:nFolds] = lapply(1:nFolds,citrus.removeFeatures,foldFeatures=foldFeatures,nonOverlappingFeatures=nof)
#    leftoutFeatures = lapply(1:nFolds,citrus.removeFeatures,foldFeatures=leftoutFeatures,nonOverlappingFeatures=nof)
#    conditionResults = list(foldLargeEnoughClusters=foldLargeEnoughClusters,foldFeatures=foldFeatures,leftoutFeatures=leftoutFeatures,folds=folds)
#  } else {
#    conditionResults = list(foldLargeEnoughClusters=foldLargeEnoughClusters,foldFeatures=foldFeatures,folds=folds)  
#  }
#  save(conditionResults,file=file.path(outputDir,paste("citrus.Features.",conditionName,".rDat",sep="")))
#  return(conditionResults)
#}

#citrus.buildFoldFeatures = function(index,featureTypes=c("abundances"),folds,citrus.dataArray,foldsClusterAssignments,foldLargeEnoughClusters,conditions,calculateLeaveoutData=F,...){
#  if ((length(folds[[index]])==1) && (folds[[index]] %in% c("all","mappingResult"))){
#    foldsFileIds=as.vector(citrus.dataArray$fileIds[,conditions])
#  } else if (calculateLeaveoutData){
#    foldsFileIds=as.vector(citrus.dataArray$fileIds[folds[[index]],conditions])
#  } else {
#    foldsFileIds=as.vector(citrus.dataArray$fileIds[-folds[[index]],conditions])
#  }
#  return(citrus.buildClusterFeatures(clusterAssignments=foldsClusterAssignments[[index]],featureTypes,largeEnoughClusters=foldLargeEnoughClusters[[index]],foldsFileIds=foldsFileIds,conditions=conditions,citrus.dataArray=citrus.dataArray,...))
#}

#citrus.buildClusterFeatures = function(clusterAssignments,featureTypes,largeEnoughClusters,foldsFileIds,conditions,citrus.dataArray,...){
#  features = list()
#  for (featureType in sort(featureTypes)){
#     features[[featureType]] = do.call(paste("citrus.calculateFeature",featureType,sep="."),args=list(foldsFileIds=foldsFileIds,clusterIds=largeEnoughClusters,clusterAssignments=clusterAssignments,data=citrus.dataArray$data[(citrus.dataArray$data[,"fileId"]%in%foldsFileIds),],conditions=conditions,citrus.dataArray=citrus.dataArray,preCalcFeatures=features,...=...))
#  }
#  for (featureType in sort(featureTypes)){
# ASSUME THAT WANT FEATURE DIFFERENCES. MAY BE BAD ASSUMPTION
#    if ((length(conditions)==2) && (featureType!="emDists")){
#      cat("Two conditions found. Assuming differential features of interest.\n")
#      cat(paste("Caluclating difference in ",featureType," between ",conditions[2]," & ",conditions[1],".\n",sep=""))
#      fns2 = citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[2]]]
#      fns1 = citrus.dataArray$fileNames[citrus.dataArray$fileIds[,conditions[1]]]
#      features[[paste(featureType,"difference")]] = features[[featureType]][rownames(features[[featureType]]) %in% fns2,] - features[[featureType]][rownames(features[[featureType]]) %in% fns1,]
#      colnames(features[[paste(featureType,"difference")]]) = paste(colnames(features[[paste(featureType,"difference")]]),"difference")
#      features[[featureType]]=NULL
#    }
#    
#  }
#  return(do.call("cbind",features))
#}




#citrus.calculateFeature.emDists = function(foldsFileIds,clusterIds,clusterAssignments,data,foldFileNames,conditions,citrus.dataArray,preCalcFeatures,...){
#  library("emdist")
#  addtlArgs = list(...)
#  if (!("emdColumns" %in% names(addtlArgs))){
#    stop("emdColumns argument must be specified to compute cluster emDists.")
#  }
#  if (length(conditions)!=2){
#    warning("Only know how to calculate EMD Features for two conditions.")
#    return(NULL)
#  }
#  if ("abundances" %in% names(preCalcFeatures)){
#    df = preCalcFeatures[["abundances"]]
#  } else {
#    df = citrus.calculateFeature.abundances(foldsFileIds,clusterIds,clusterAssignments,data,citrus.dataArray,...)
#  }
#  completeClusterIds = clusterIds[apply(df>0,2,all)]

#  referenceFileIds = foldsFileIds[foldsFileIds %in% citrus.dataArray$fileIds[,conditions[1]]]
#  targetFileIds = foldsFileIds[foldsFileIds %in% citrus.dataArray$fileIds[,conditions[2]]]
#  if ("mc.cores" %in% names(addtlArgs)){
#    features = do.call("rbind",mclapply(1:length(referenceFileIds),citrus.calculateFileClustersEMDist,clusterIds=completeClusterIds,clusterAssignments=clusterAssignments,referenceFileIds=referenceFileIds,targetFileIds=targetFileIds,data=data,emdColumns=addtlArgs[["emdColumns"]],mc.cores=addtlArgs[["mc.cores"]]))
#  } else {
#    features = t(sapply(1:length(referenceFileIds),citrus.calculateFileClustersEMDist,clusterIds=completeClusterIds,clusterAssignments=clusterAssignments,referenceFileIds=referenceFileIds,targetFileIds=targetFileIds,data=data,emdColumns=addtlArgs[["emdColumns"]]))
#  }
#  rownames(features) = citrus.dataArray$fileNames[targetFileIds]
#  return(features)
#}

#citrus.calculateFileClustersEMDist = function(sampleIndex,clusterIds,clusterAssignments,referenceFileIds,targetFileIds,data,emdColumns){
#  res = sapply(clusterIds,citrus.calculateFileClusterEMDist,clusterAssignments=clusterAssignments,referenceFileId=referenceFileIds[sampleIndex],targetFileId=targetFileIds[sampleIndex],data=data,emdColumns=emdColumns)
#  res2 = as.vector(res) 
#  names(res2) = paste(paste(paste("cluster",rep(clusterIds,each=nrow(res))),rownames(res)),"emDist")
#  return(res2)  
#}

#citrus.calculateFileClusterEMDist = function(clusterId,clusterAssignments,referenceFileId,targetFileId,data,emdColumns){
#  clusterData = data[clusterAssignments[[clusterId]],]
#  referenceData = clusterData[clusterData[,"fileId"]==referenceFileId,]
#  targetData = clusterData[clusterData[,"fileId"]==targetFileId,]
#  res=sapply(emdColumns,citrus.calculateFileClusterParameterEMDist,referenceData,targetData)
#  if (is.numeric(emdColumns)){
#    names(res)=colnames(data)[emdColumns]  
#  } else {
#    names(res)=emdColumns 
#  }
#  return(res)
#}

#citrus.calculateFileClusterParameterEMDist = function(emdColumn,referenceData,targetData){
#cat("IMPLEMENT MINIMUM CLUSTER PERCENTAGE CHECK\n");
#  if ((length(referenceData)==0)||(length(targetData)==0)){
#    return(0);
#  } else if (is.null(dim(referenceData))||is.null(dim(targetData))){
#    return(0);
#  } else if ((nrow(referenceData)<25)||(nrow(targetData)<25)){
#    return(0);
#  } 
#  h = hist(c(referenceData[,emdColumn],targetData[,emdColumn]),breaks=50,plot=F)
#  dist = emdw(A=h$mids,wA=hist(referenceData[,emdColumn],plot=F,breaks=h$breaks)$density,B=h$mids,wB=hist(targetData[,emdColumn],plot=F,breaks=h$breaks)$density)
#  if ((mean(referenceData[,emdColumn])-mean(targetData[,emdColumn]))>0){
#    dist = -dist
#  }
#  return(dist)
#}

#citrus.calculateFeature.abundances = function(foldsFileIds,clusterIds,clusterAssignments,data,citrus.dataArray,...){
#  addtlArgs = list(...)
#  if ("mc.cores" %in% names(addtlArgs)){
#    features = do.call("rbind",mclapply(foldsFileIds,citrus.calculateFileClustersAbundances,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))
#  } else {
#    features = t(sapply(foldsFileIds,citrus.calculateFileClustersAbundances,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))
#  }
#  #features = t(sapply(foldsFileIds,citrus.calculateFileClustersAbundances,clusterIds=clusterIds,clusterAssignments=clusterAssignments,data=data,...))
#  rownames(features) = citrus.dataArray$fileNames[foldsFileIds]
#  return(features)
#}

#citrus.calculateFileClustersAbundances = function(fileId,clusterIds,clusterAssignments,data,...){
#  fileIds=data[,"fileId"]
#  res = sapply(clusterIds,citrus.calculateFileClusterAbundance,clusterAssignments=clusterAssignments,fileId=fileId,fileIds=fileIds)
#  names(res) = paste(paste("cluster",clusterIds),"abundance")
#  return(res)
#}

#citrus.calculateFileClusterAbundance = function(clusterId,clusterAssignments,fileId,fileIds){
#sum(which(fileIds==fileId) %in% clusterAssignments[[clusterId]])/sum((fileIds==fileId))
#}


