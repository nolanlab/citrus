citrus.plotTypeErrorRate = function(modelType,outputDir,regularizationThresholds,thresholdCVRates,foldModels,cvMinima,family){  
    nAllFolds = length(foldModels[[modelType]])
    modelOutputDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    dir.create(modelOutputDir)
    pdf(file.path(modelOutputDir,"ModelErrorRate.pdf"),width=6,height=6)
    thresholds=regularizationThresholds[[modelType]]
    errorRates=thresholdCVRates[[modelType]][,"cvm"]
    ylim=c(0,1)
    ylab="Model Cross Validation Error Rate"
    if (modelType=="glmnet"){
      thresholds = log(thresholds)
      xlab="log(Regularization Threshold)"
      nonzeroCounts = foldModels[[modelType]][[nAllFolds]]$df
      if (family=="survival"){
        ylim=range(errorRates,na.rm=T)
        ylab="Model Cross Validation Partial Likelihood"
      }
    } else if (modelType=="pamr"){
      xlab="Regularization Threshold"
      nonzeroCounts = foldModels[[modelType]][[nAllFolds]]$nonzero
    }
    plot(errorRates,type='o',pch=20,col="red",main="Number of model features\n",axes=F,xlab=xlab,ylim=ylim,ylab=ylab)
    #Plot SEM
    for (i in 1:length(thresholds)){
      lines(c(i,i),c(errorRates[i]+thresholdCVRates[[modelType]][,"cvsd"][i],errorRates[i]-thresholdCVRates[[modelType]][,"cvsd"][i]),col="red",lty=3)
    }
    grid()
    axis(1,at=1:length(errorRates),labels=sapply(thresholds,citrus.formatDecimal))
    if (family=="survival"){
      axis(2)
    } else {
      axis(2,at=c(0,.25,.5,.75,1),labels=c(0,25,50,75,100))  
    }
    
    axis(3,labels=nonzeroCounts,at=1:length(errorRates))
    legendLabels = c("Cross Validation Error Rate")
    legendColors = c("red")
    legendPchs = c(20)
    legendLty = c(1)
    legendPtCex=c(1)
    if (!is.null(thresholdCVRates[[modelType]]$fdr)){
      lines(thresholdCVRates[[modelType]]$fdr,type='o',pch=1,cex=1.5,col="blue")
      legendLabels = c(legendLabels,"Feature False Discovery Rate")
      legendColors = c(legendColors,"blue")
      legendPchs = c(legendPchs,1)
      legendLty = c(legendLty,1)
      legendPtCex=c(legendPtCex,1)
    }
    
    cv.min = cvMinima[[modelType]]$cv.min.index
    cv.1se = cvMinima[[modelType]]$cv.1se.index
    if (!is.null(cv.min)){
      points(c(cv.min,cv.min),y=c(errorRates[cv.min],errorRates[cv.min]),col="green",pch=20,cex=2)  
    }
    if (!is.null(cv.1se)){
      points(c(cv.1se,cv.1se),y=c(errorRates[cv.1se],errorRates[cv.1se]),col="orange",pch=9,cex=2)  
    }
        
    legendLabels = c(legendLabels,"cv.min","cv.1se")
    legendColors = c(legendColors,"green","orange")
    legendPchs = c(legendPchs,20,9)
    legendLty = c(legendLty,0,0)
    legendPtCex=c(legendPtCex,2,1.5)  
       
    if ( "cv.fdr.constrained" %in% names(cvMinima[[modelType]])) {
      
      cv.fdr.constrained = cvMinima[[modelType]]$cv.fdr.constrained.index
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="yellow",pch=17,cex=1.5)
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="black",pch=2,cex=1.5)
      legendLabels = c(legendLabels,"cv.fdr.constrained")
      legendColors = c(legendColors,"yellow")
      legendPchs = c(legendPchs,17)
      legendLty = c(legendLty,0)
      legendPtCex=c(legendPtCex,1.5)  
      
    }
    legend(x="topleft",legendLabels,col=legendColors,pch=legendPchs,lty=legendLty,pt.cex=legendPtCex,cex=.8,bg="white")
    dev.off()   
}

citrus.plotDifferentialFeatures = function(modelType,differentialFeatures,features,outputDir,labels,family,foldModels,cvMinima,regularizationThresholds){
  do.call(paste("citrus.plotModelDifferentialFeatures",family,sep="."),args=list(modelType=modelType,differentialFeatures=differentialFeatures,features=features,outputDir=outputDir,labels=labels,foldModels=foldModels,cvMinima=cvMinima,regularizationThresholds=regularizationThresholds))
}

citrus.plotModelDifferentialFeatures.survival = function(modelType,differentialFeatures,features,outputDir,labels,foldModels,cvMinima,regularizationThresholds,...){
  s = Surv(time=labels[,"time"],event=labels[,"event"])
  for (cvPoint in names(differentialFeatures[[modelType]])){
    modelTypeDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    nonzeroFeatureNames = differentialFeatures[[modelType]][[cvPoint]][["features"]]
    pchs = rep(13,nrow(s))
    pchs[as.logical(s[,2])]=20
    pdf(file.path(modelTypeDir,paste("features_",cvPoint,".pdf",sep="")))
    for (nonzeroFeatureName in nonzeroFeatureNames){
      f = features[,nonzeroFeatureName]
      cutoff = median(f)
      plot(x=f,y=s[,1],xlab=nonzeroFeatureName,ylab="time",pch=pchs,col="red",cex=2,main="Time vs Feature")   
      #lines(c(cutoff,cutoff),c(min(s[,1]-10),max(s[,1]+10)),col="blue",lty=3)
      legend("topright",legend=c("censored","uncensored","median"),pch=c(13,20,45),col=c("red","red","blue"))
    }
    dev.off();
    pdf(file.path(modelTypeDir,paste("survivalCurves_singleFeatures_",cvPoint,".pdf",sep="")),height=6,width=6)
    for (nonzeroFeatureName in nonzeroFeatureNames){
      f = features[,nonzeroFeatureName]
      cutoff = median(f)
      sf=survfit(s~group,data=data.frame(f,group=as.numeric(f<cutoff)))
      plot(sf,xlab="Time",ylab="Percent Survival",main=paste("Survival stratified on",nonzeroFeatureName))
      sdf = survdiff(s~group,data=data.frame(f,group=as.numeric(f<cutoff)))
      p.val = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      legend("topright",legend=c(paste("P-Value:",substr(p.val,1,5))),cex=1)
    }
    dev.off()
    finalModel=foldModels[[modelType]][[length(foldModels[[modelType]])]]
    pdf(file.path(modelTypeDir,paste("survivalCurves_allFeatures_",cvPoint,".pdf",sep="")))
    f=citrus.predict.survival(finalModel,features,s=cvMinima[[modelType]][[cvPoint]])
    cutoff = median(f)
    sf=survfit(s~group,data=data.frame(f,group=as.numeric(f<cutoff)))
    plot(sf,xlab="Time",ylab="Percent Survival",main=paste("Survival stratified on",cvPoint,"model"))
    sdf = survdiff(s~group,data=data.frame(f,group=as.numeric(f<cutoff)))
    legend("topright",legend=c(1 - pchisq(sdf$chisq, length(sdf$n) - 1)),cex=.7)
    dev.off()
  }
}

citrus.plotModelDifferentialFeatures.classification = function(modelType,differentialFeatures,features,outputDir,labels,...){
  for (cvPoint in names(differentialFeatures[[modelType]])){
    modelTypeDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    nonzeroFeatureNames = differentialFeatures[[modelType]][[cvPoint]][["features"]]
    for (nonzeroFeatureName in nonzeroFeatureNames){
      df = data.frame(value=features[,nonzeroFeatureName],labels=as.factor(labels),featureName=nonzeroFeatureName)
      if (which(nonzeroFeatureNames==nonzeroFeatureName)==1){
        combinedDf = df;
      } else {
        combinedDf = rbind(combinedDf,df)
      }  
    }
    nrow=ceiling(length(nonzeroFeatureNames)/4)
    ncol=4
    if (length(nonzeroFeatureNames)<4){
      ncol=length(nonzeroFeatureNames)
    }
    
    pdf(file.path(modelTypeDir,paste("features-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=(ncol*4),height=(nrow*1.5))
    p <- ggplot(combinedDf[,], aes(labels, value)) 
    p = p + facet_wrap(~featureName,ncol=4) + geom_boxplot(outlier.colour=rgb(0,0,0,0),colour=rgb(0,0,0,.3)) + geom_point(aes(color=labels),alpha=I(0.25),shape=19,size=I(2)) + scale_colour_manual(values = c("red","blue")) + coord_flip() +  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none") 
    print(p)
    dev.off()
  }
  
}

scaleToRange =function(x,scale){
  xrange = max(x)-min(x)
  x = x/xrange
  x = x-min(x)
  scaleRange = max(scale)-min(scale)
  x=x*scaleRange
  x = x+min(scale)
  return(x)
}

citrus.overlapDensityPlot = function(clusterDataList,backgroundData){
  combined = data.frame(check.names=F,check.rows=F)
  for (clusterName in names(clusterDataList)){
    combined=rbind(combined,data.frame(value=as.vector(clusterDataList[[clusterName]]),marker=as.vector(sapply(colnames(clusterDataList[[clusterName]]),rep,nrow(clusterDataList[[clusterName]]))),clusterId=clusterName,dplot="cluster",check.names=F,check.rows=F))
    combined=rbind(combined,data.frame(value=as.vector(backgroundData),marker=as.vector(sapply(colnames(clusterDataList[[clusterName]]),rep,nrow(backgroundData))),clusterId=clusterName,dplot="background",check.names=F,check.rows=F))
  }
  p = ggplot(combined) + geom_density(aes(x=value,fill=dplot,colour=dplot),alpha=.6) + facet_grid(clusterId~marker,scales="free")
  print(p)
}

citrus.plotClusters = function(modelType,differentialFeatures,outputDir,clusterChildren,citrus.dataArray,conditions,clusterCols,clusterColLabels=NULL){
  data = citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions],]
  clusterChildren = clusterChildren[[length(clusterChildren)]]
  if (!is.null(clusterColLabels)){
    if (length(clusterColLabels)!=length(clusterCols)){
      stop("clusterColsLabels length must equal length of clusterCols");
    }
  }
  for (cvPoint in names(differentialFeatures[[modelType]])){
    nonzeroClusters = as.numeric(differentialFeatures[[modelType]][[cvPoint]][["clusters"]])
    pdf(file=file.path(outputDir,paste(modelType,"_results/clusters-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=(2.2*length(clusterCols)+2),height=(2.2*length(nonzeroClusters)))
    clusterDataList=list();
    for (nonzeroCluster in sort(nonzeroClusters)){
      if (nrow(data[clusterChildren[[nonzeroCluster]],])>5000){
        clusterDataList[[as.character(nonzeroCluster)]]=data[clusterChildren[[nonzeroCluster]],clusterCols][sample(1:nrow(data[clusterChildren[[nonzeroCluster]],]),1000),]
      } else {
        clusterDataList[[as.character(nonzeroCluster)]]=data[clusterChildren[[nonzeroCluster]],clusterCols]
      }
      if (!is.null(clusterColLabels)){
        colnames(clusterDataList[[as.character(nonzeroCluster)]])=clusterColLabels	
      }
    }
    if (nrow(data)>2500){
      bgData = data[sample(1:nrow(data),2500),clusterCols]
    } else {
      bgData = data[,clusterCols]
    }
    
    citrus.overlapDensityPlot(clusterDataList=clusterDataList,backgroundData=bgData)
    dev.off()
  }
}

########################################
# Hierarchical Clustering plots
########################################

.getClusterMedians = function(clusterId,clusterAssignments,clusterCols,data){
  apply(data[clusterAssignments[[clusterId]],clusterCols],2,median)
}

.decimalFormat = function(x){
  sprintf("%.2f",x)
}

.scaleToOne = function(x){
  x = x-min(x)
  x = x/max(x)
  return(x)
}

.petalVertex <- function(coords, params){
  scale = params("vertex", "scale")
  weights = params("vertex", "weights")
  sapply(1:nrow(coords),.petalVertexWrapper,coordinates=coords,scale=scale,weights=weights)
}

.petalVertexWrapper = function(x,coordinates,scale,weights){
  .petalPlot(xpos=coordinates[x,1],coordinates[x,2],d=weights[x,],scale=scale)
}

.petalPlot = function(xpos,ypos,d,scale=1,segCol=NULL,labels=NULL){
  res=360
  nSegments = length(d)
  if (is.null(segCol)){
    segCol=rainbow(nSegments,alpha=.7)
    borCol=rainbow(nSegments)
  }
  angles = (360/nSegments)*c(0:nSegments)*(pi/180)
  series = seq(from=1,to=360,length.out=(res))
  points = series*(pi/180)
  x = cos(points)
  y = sin(points)
  segPoints = floor(seq(from=1,to=360,length.out=(nSegments+1)))
  for (i in 1:nSegments){
    seg.x = c(xpos,xpos+x[segPoints[i]:segPoints[i+1]]*d[i]*scale)
    seg.y = c(ypos,ypos+y[segPoints[i]:segPoints[i+1]]*d[i]*scale)
    polygon(seg.x,seg.y,col=segCol[i],border=borCol[i])
    lines(c(xpos,xpos+x[segPoints[i]]*scale),c(ypos,ypos+y[segPoints[i]]*scale),col=rgb(1,1,1,.4))
  }
  lines((x*scale)+xpos,(y*scale)+ypos,col=rgb(1,1,1,.4))
  lines((x*.5*scale)+xpos,(y*.5*scale)+ypos,col=rgb(1,1,1,.4))
  
  if (!is.null(labels)){
    labelPoints = segPoints[-length(segPoints)]+segPoints[2]/2
    text(x=c(x[labelPoints]*scale+xpos),y=c(y[labelPoints]*scale+ypos),labels,col="white")
  }
  
}

.getClusterFeatureMatrix = function(featureVector){
  df = do.call("rbind",strsplit(gsub(pattern="(cluster [0-9]+) ",replacement="\\1~",featureVector),"~"))
  return(cbind(cluster=do.call("rbind",strsplit(df[,1]," "))[,2],feature=df[,2]))
}
  

citrus.createHierarchyGraph = function(largeEnoughClusters,mergeOrder,clusterAssignments){
  cmat = matrix(0,ncol=length(largeEnoughClusters),nrow=length(largeEnoughClusters),dimnames=list(largeEnoughClusters,largeEnoughClusters))
  for (cluster in largeEnoughClusters){
    children = mergeOrder[cluster,]  
    for (child in children){
      if (as.character(child) %in% rownames(cmat)){
        cmat[as.character(cluster),as.character(child)] = 1
      }
    }
  }
  g = graph.adjacency(adjmatrix=cmat,mode="directed",add.colnames='label')
  clusterSizes = do.call("rbind",lapply(clusterAssignments,length))
  sizes = sapply(log(clusterSizes[largeEnoughClusters]),findInterval,vec=hist(log(clusterSizes[largeEnoughClusters]),breaks=10,plot=F)$breaks)
  sizes = sizes+6
  g = set.vertex.attribute(g,"size",value=sizes)
  return(g)
}

citrus.plotHierarchicalClusterMedians = function(outputFile,clusterMedians,graph,layout){
  pdf(file=outputFile,width=15,height=15,bg="black")
  for (target in 1:ncol(clusterMedians)){
    ct = seq(from=(min(clusterMedians[,target])-0.01),to=(max(clusterMedians[,target])+0.01),length.out=20)
    cols = topo.colors(20)[sapply(clusterMedians[,target],findInterval,vec=ct)]
    par(col.main="white")  
    plot.igraph(graph,layout=layout,vertex.color=cols,main=colnames(clusterMedians)[target],edge.color="white",vertex.label.color="white",edge.arrow.size=.2,vertex.frame.color=rgb(1,1,1,.5),vertex.label.cex=.7,vertex.label.family="Helvetica")
    
    # Legend
    legend_image <- as.raster(matrix(rev(topo.colors(20)), ncol=1))
    rasterImage(legend_image, 1.1, -.5, 1.15,.5)
    text(x=1.15, y = seq(-.5,.5,l=5), labels = .decimalFormat(ct[c(1,floor((length(ct)/4)*1:4))]) ,pos=4,col="white")
  }
  dev.off()  
}


citrus.plotHierarchicalClusterFeatureGroups = function(outputFile,featureClusterMatrix,largeEnoughClusters,graph,layout,petalPlots=F,clusterMedians=NULL){
  pdf(file=outputFile,width=15,height=15,bg="black")
  for (feature in unique(featureClusterMatrix[,"feature"])){
    fGroup = list();
    featureClusters = as.numeric(featureClusterMatrix[featureClusterMatrix[,"feature"]==feature,"cluster"])
    featureElements = match(featureClusters,largeEnoughClusters)
    subgraph = induced.subgraph(graph,featureElements)
    groupAssignments = clusters(subgraph)$membership
    for (groupId in unique(groupAssignments)){
      fGroup[[groupId]]=match(get.vertex.attribute(subgraph,"label")[groupAssignments==groupId],get.vertex.attribute(graph,"label"))
    }
    par(col.main="white")
    if (petalPlots){
      if (is.null(clusterMedians)){
        stop("clusterMedians argument must be supplied to plot petals")
      }
      add.vertex.shape("petal", clip=vertex.shapes("circle")$clip,plot=.petalVertex, parameters=list(vertex.scale=.04,vertex.weights=apply(clusterMedians,2,.scaleToOne)))
      plot.igraph(graph,layout=layout,mark.groups=fGroup,mark.expand=5,main=feature,edge.color="white",vertex.label.color="white",edge.arrow.size=.2,vertex.frame.color=rgb(1,1,1,.5),vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=rgb(0,0,.5,.3),mark.col=topo.colors(length(fGroup),alpha=.3),vertex.shape="petal")
      .petalPlot(xpos=1,ypos=-1,d=rep(1,ncol(clusterMedians)),scale=.2,labels=colnames(clusterMedians))
    } else {
      plot.igraph(graph,layout=layout,mark.groups=fGroup,mark.expand=5,main=feature,edge.color="white",vertex.label.color="white",edge.arrow.size=.2,vertex.frame.color=rgb(1,1,1,.5),vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=rgb(0,0,.5,.5),mark.col=topo.colors(length(fGroup),alpha=.8))    
    }
  }
  dev.off()
}




# Plot
citrus.plotRegressionResults = function(outputDir,citrus.preclusterResult,citrus.featureObject,citrus.regressionResult,modelTypes,family,labels,plotTypes=c("errorRate","stratifyingFeatures","stratifyingClusters","clusterGraph"),...){
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
    
    if ("errorRate" %in% plotTypes){
      cat("Plotting Error Rate\n")
      sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds,thresholdCVRates=citrus.regressionResult[[conditionName]]$thresholdCVRates,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,family=family)
    }
    
    if ("stratifyingFeatures" %in% plotTypes){
      cat("Plotting Stratifying Features\n")
      lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,features=citrus.featureObject[[conditionName]]$foldFeatures[[nAllFolds]],outputDir=conditionOutputDir,labels=labels,family=family,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds)
    }
    
    if ("stratifyingClusters" %in% plotTypes){
      cat("Plotting Stratifying Clusters\n")
      lapply(modelTypes,citrus.plotClusters,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,outputDir=conditionOutputDir,clusterChildren=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments,citrus.dataArray=citrus.preclusterResult[[conditionName]]$citrus.dataArray,conditions=citrus.preclusterResult[[conditionName]]$conditions,clusterCols=citrus.preclusterResult[[conditionName]]$clusterColumns,clusterColLabels=clusterColLabels)
    }
    
    
    # GO BACK AND MAKE THIS PARALLEL CALLS....
    if ("clusterGraph" %in% plotTypes){
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
  
}
