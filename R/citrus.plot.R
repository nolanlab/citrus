citrus.plotTypeErrorRate = function(modelType,outputDir,regularizationThresholds,thresholdCVRates,foldModels,cvMinima,family){  
    if (modelType=="sam"){
      return()
    }
    nAllFolds = length(foldModels[[modelType]])
    modelOutputDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
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
  # Passing sam models makes things go crazy so this is a hack to avoid problems debugging. 
  if (family=="survival"){
    do.call(paste("citrus.plotModelDifferentialFeatures",family,sep="."),args=list(modelType=modelType,differentialFeatures=differentialFeatures,features=features,outputDir=outputDir,labels=labels,foldModels=foldModels,cvMinima=cvMinima,regularizationThresholds=regularizationThresholds))    
  } else {
    do.call(paste("citrus.plotModelDifferentialFeatures",family,sep="."),args=list(modelType=modelType,differentialFeatures=differentialFeatures,features=features,outputDir=outputDir,labels=labels,cvMinima=cvMinima,regularizationThresholds=regularizationThresholds))
  }
}

citrus.plotModelDifferentialFeatures.survival = function(modelType,differentialFeatures,features,outputDir,labels,foldModels,cvMinima,regularizationThresholds,...){
  s = Surv(time=labels[,"time"],event=labels[,"event"])
  for (cvPoint in names(differentialFeatures[[modelType]])){
    modelTypeDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    nonzeroFeatureNames = differentialFeatures[[modelType]][[cvPoint]][["features"]]
    
    # Write features to file for easy parsing
    write.table(features[,nonzeroFeatureNames],file=file.path(modelTypeDir,paste("features_",cvPoint,".csv",sep="")),quote=F,sep=",")
    
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

    # Write features to file for easy parsing
    write.table(features[,nonzeroFeatureNames],file=file.path(modelTypeDir,paste("features_",cvPoint,".csv",sep="")),quote=F,sep=",")

    melted = melt(data.frame(features[,nonzeroFeatureNames,drop=F],labels=labels,check.names=F),id.vars="labels")
    
    nrow=ceiling(length(nonzeroFeatureNames)/4)
    ncol=4
    if (length(nonzeroFeatureNames)<4){
      ncol=length(nonzeroFeatureNames)
    }
    
    pdf(file.path(modelTypeDir,paste("features-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=(ncol*4),height=(nrow*1.5))
    p <- ggplot(melted, aes(x=factor(labels), y=value)) 
    p = p + facet_wrap(~variable,ncol=4) + geom_boxplot(outlier.colour=rgb(0,0,0,0),colour=rgb(0,0,0,.3)) + geom_point(aes(color=factor(labels)),alpha=I(0.25),shape=19,size=I(2)) + coord_flip() +  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none") 
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
    combined=rbind(combined,cbind(melt(clusterDataList[[clusterName]],varnames=c("row","marker")),clusterId=clusterName,src="Cluster"))
  }
  p = ggplot(data=combined,aes(x=value, y = ..scaled..,fill=src)) + geom_density() + facet_grid(clusterId~marker,scales="free")+geom_density(data=cbind(melt(backgroundData,varnames=c("row","marker")),src="Background"))+theme_bw()+scale_fill_manual(values = c("Background" = rgb(.3,.3,1,.2), "Cluster" = rgb(1,.3,.3,.5)))
  p = p+theme(legend.position="bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_blank())
  p = p+labs(fill="Distribution:")
  print(p)
}

citrus.plotModelClusters = function(modelType,differentialFeatures,outputDir,clusterAssignments,citrus.dataArray,conditions,clusterCols){
  for (cvPoint in names(differentialFeatures[[modelType]])){
    clusterIds = as.numeric(differentialFeatures[[modelType]][[cvPoint]][["clusters"]])
    outputFile = file.path(outputDir,paste(modelType,"_results/clusters-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep=""))
    citrus.plotClusters(clusterIds,clusterAssignments=clusterAssignments[[length(clusterAssignments)]],citrus.dataArray,conditions,clusterCols,outputFile=outputFile)
  }
}

citrus.plotClusters = function(clusterIds,clusterAssignments,citrus.dataArray,conditions,clusterCols,outputFile=NULL){
  data = citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions],]
  if (!is.null(outputFile)){
    pdf(file=outputFile,width=(2.2*length(clusterCols)+2),height=(2*length(clusterIds)))  
  }
  clusterDataList = list();
  for (clusterId in sort(clusterIds)){
    if (length(clusterAssignments[[clusterId]])>2500){
      clusterDataList[[as.character(clusterId)]]=data[clusterAssignments[[clusterId]],clusterCols][sample(1:length(clusterAssignments[[clusterId]]),2500),]
    } else {
      clusterDataList[[as.character(clusterId)]]=data[clusterAssignments[[clusterId]],clusterCols]
    }
    
    colLabels = citrus.dataArray$fileChannelNames[[conditions[1]]][[1]]
    reagentNames = citrus.dataArray$fileReagentNames[[conditions[1]]][[1]]
    displayNames = colLabels
    displayNames[nchar(reagentNames)>1] = reagentNames[nchar(reagentNames)>1]
    if (is.numeric(clusterCols)){
      colnames(clusterDataList[[as.character(clusterId)]])=displayNames[clusterCols]  
    } else {
      colnames(clusterDataList[[as.character(clusterId)]])=displayNames[colLabels%in%clusterCols]
    }
  }
  if (nrow(data)>2500){
    bgData = data[sample(1:nrow(data),2500),clusterCols]
  } else {
    bgData = data[,clusterCols]
  }
  colnames(bgData) = colnames(clusterDataList[[1]])
  citrus.overlapDensityPlot(clusterDataList=clusterDataList,backgroundData=bgData)
  if (!is.null(outputFile)){
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



.getClusterFeatureMatrix = function(featureVector){
  df = do.call("rbind",strsplit(gsub(pattern="(cluster [0-9]+) ",replacement="\\1~",featureVector),"~"))
  return(cbind(cluster=do.call("rbind",strsplit(df[,1]," "))[,2],feature=df[,2]))
}
  

citrus.createHierarchyGraph = function(largeEnoughClusters,mergeOrder,clusterAssignments,minVertexSize=6){
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
  sizes = sizes+minVertexSize
  g = set.vertex.attribute(g,"size",value=sizes)
  return(g)
}

citrus.plotHierarchicalClusterMedians = function(outputFile,clusterMedians,graph,layout,theme="black",plotSize=15,singlePDF=F,ncol=3,scale=1,plotClusterIDs=T){
  if (theme=="black"){
    bg="black"
    stroke="white"
    strokea=rgb(1,1,1,.5)
  } else if (theme=="white"){
    bg="white"
    stroke="black"
    strokea=rgb(0,0,0,.5)
  } else {
    stop("Unrecognized theme option. Choices are 'white' or 'black'")
  }
  if (plotClusterIDs){
    vc="white"
  } else {
    vc=rgb(0,0,0,0)
  }
  if (singlePDF){
    expand=3*scale
    nrow = ceiling(ncol(clusterMedians)/ncol)
    pdf(file=outputFile,width=ncol*expand,height=nrow*expand,bg=bg)
    par(mfrow=c(nrow,ncol),oma=rep(0.8,4),mar=c(.1,0,.8,2.8))
  } else {
    pdf(file=outputFile,width=plotSize,height=plotSize,bg=bg)  
    
  }
  
  for (target in 1:ncol(clusterMedians)){
    ct = seq(from=(min(clusterMedians[,target])-0.01),to=(max(clusterMedians[,target])+0.01),length.out=20)
    cols = .graphColorPalette(20)[sapply(clusterMedians[,target],findInterval,vec=ct)]
    par(col.main=stroke)  
    plot.igraph(graph,layout=layout,vertex.color=cols,main=colnames(clusterMedians)[target],edge.color=stroke,vertex.label.color=vc,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica")
        
    # Legend
    legend_image <- as.raster(matrix(rev(.graphColorPalette(20)), ncol=1))
    rasterImage(legend_image, 1.1, -.5, 1.15,.5)
    text(x=1.15, y = seq(-.5,.5,l=5), labels = .decimalFormat(ct[c(1,floor((length(ct)/4)*1:4))]) ,pos=4,col=stroke)
  }
  dev.off()  
}


citrus.plotHierarchicalClusterFeatureGroups = function(outputFile,featureClusterMatrix,graph,layout,clusterMedians=NULL,featureClusterCols=NULL,theme="black",encircle=T,plotSize=15,plotClusterIDs=T){
  
  if (!is.null(featureClusterCols)&&(is.null(names(featureClusterCols)))){
    stop("featureClusterCols argument must be vector with elements having names of vertices to be colored.")
  }
  
  if (theme=="black"){
    bg="black"
    stroke="white"
    strokea=rgb(1,1,1,.5)
  } else if (theme=="white"){
    bg="white"
    stroke="black"
    strokea=rgb(0,0,0,.5)
  } else {
    stop("Unrecognized theme option. Choices are 'white' or 'black'")
  }
  pdf(file=outputFile,width=plotSize,height=plotSize,bg=bg)
  for (feature in unique(featureClusterMatrix[,"feature"])){
    fGroup = list();
    featureClusters = featureClusterMatrix[featureClusterMatrix[,"feature"]==feature,"cluster"]
    featureElements = match(featureClusters,as.numeric(get.vertex.attribute(graph,"label",V(graph))))
    subgraph = induced.subgraph(graph,featureElements)
    groupAssignments = clusters(subgraph)$membership
    for (groupId in unique(groupAssignments)){
      fGroup[[groupId]]=match(get.vertex.attribute(subgraph,"label")[groupAssignments==groupId],get.vertex.attribute(graph,"label"))
    }
    par(col.main=stroke)

    vertexFont=rep(1,length(V(graph)))
    vertexFont[get.vertex.attribute(graph,"label")%in%featureClusters]=2
    if (is.null(featureClusterCols)){
      vertexColor=rep(rgb(0,0,.5,.5),length(V(graph)))
      vertexColor[get.vertex.attribute(graph,"label")%in%featureClusters]=rgb(0.5,0,0,.7)
    } else {
      vertexColor=rep(rgb(0,0,.5,.5),length(V(graph)))
      cp = rgb(1,0,0,seq(0,1,by=.05))
      ct = seq(from=(min(featureClusterCols)-0.01),to=(max(featureClusterCols)+0.01),length.out=20)
      vertexColor[ match(names(featureClusterCols),get.vertex.attribute(graph,"label")) ] = cp[sapply(featureClusterCols,findInterval,vec=ct)]
    }
    
    vertexLabelColor="white"
    if (!plotClusterIDs){
      vertexLabelColor = rgb(0,0,0,0)
    }
    if (encircle){
      plot.igraph(graph,layout=layout,mark.groups=fGroup,mark.expand=5,main=feature,edge.color=stroke,vertex.label.color=vertexLabelColor,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=vertexColor,mark.border=strokea,mark.col=.graphColorPalette(length(fGroup),alpha=.2),vertex.label.font=vertexFont)      
    } else {
      plot.igraph(graph,layout=layout,main=feature,edge.color=stroke,vertex.label.color=vertexLabelColor,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=vertexColor,vertex.label.font=vertexFont)     
    }
    if (!is.null(featureClusterCols)){
      legend_image <- as.raster(matrix(rev(cp), ncol=1))
      rasterImage(legend_image, 1.1, -.5, 1.15,.5)
      text(x=1.15, y = seq(-.5,.5,l=5), labels = .decimalFormat(ct[c(1,floor((length(ct)/4)*1:4))]) ,pos=4,col=stroke)  
    }
          
  }
  dev.off()
}

citrus.createPlotOutputDirectory = function(modelType,outputDir){
  modelOutputDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
  dir.create(modelOutputDir,showWarnings=F,recursive=T)
}


# Plot
citrus.plotRegressionResults = function(outputDir,citrus.preclusterResult,conditionFeatureList,citrus.regressionResult,modelTypes,family,labels,plotTypes=c("errorRate","stratifyingFeatures","stratifyingClusters","clusterGraph"),...){
  addtlArgs = list(...)
  
  theme="black"
  if ("theme" %in% names(addtlArgs)){
    theme = addtlArgs[["theme"]]
  }
  
  for (conditionName in names(citrus.regressionResult)){
    nAllFolds = length(conditionFeatureList[[conditionName]]$foldFeatures)
    # Make condition output directoy
    conditionOutputDir = file.path(outputDir,conditionName)
    sapply(modelTypes,citrus.createPlotOutputDirectory,outputDir=conditionOutputDir)
    
    
    if ("errorRate" %in% plotTypes){
      cat("Plotting Error Rate\n")
      sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=conditionOutputDir,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds,thresholdCVRates=citrus.regressionResult[[conditionName]]$thresholdCVRates,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,family=family)
    }
    
    if ("stratifyingFeatures" %in% plotTypes){
      cat("Plotting Stratifying Features\n")
      lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,features=conditionFeatureList[[conditionName]]$foldFeatures[[nAllFolds]],outputDir=conditionOutputDir,labels=labels,family=family,cvMinima=citrus.regressionResult[[conditionName]]$cvMinima,foldModels=citrus.regressionResult[[conditionName]]$foldModels,regularizationThresholds=citrus.regressionResult[[conditionName]]$regularizationThresholds)
    }
    
    if ("stratifyingClusters" %in% plotTypes){
      cat("Plotting Stratifying Clusters\n")
      lapply(modelTypes,citrus.plotModelClusters,differentialFeatures=citrus.regressionResult[[conditionName]]$differentialFeatures,outputDir=conditionOutputDir,clusterAssignments=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments,citrus.dataArray=citrus.preclusterResult[[conditionName]]$citrus.dataArray,conditions=citrus.preclusterResult[[conditionName]]$conditions,clusterCols=citrus.preclusterResult[[conditionName]]$clusterColumns)
    }
    
    
    # GO BACK AND MAKE THIS PARALLEL CALLS....
    if ("clusterGraph" %in% plotTypes){
      mcsp=0.05
      if ("mcsp" %in% names(addtlArgs)){
        mcsp=addtlArgs[["mcsp"]]
      }
      if (mcsp<0.005){
        minVertexSize=0
        plotSize=35
      } else if (mcsp<0.01){
        minVertexSize=4
        plotSize=20
      } else if (mcsp<0.05){
        minVertexSize=6
        plotSize=15
      } else {
        minVertexSize=8
        plotSize=10
      }
      
      cat("Plotting Clustering Graph\n")
      g = citrus.createHierarchyGraph(largeEnoughClusters=conditionFeatureList[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],mergeOrder=citrus.preclusterResult[[conditionName]]$foldsCluster[[nAllFolds]]$merge,clusterAssignments=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]],minVertexSize=minVertexSize)
      l = layout.reingold.tilford(g,root=length(V(g)),circular=T)
      clusterMedians = t(sapply(conditionFeatureList[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]],.getClusterMedians,clusterAssignments=citrus.preclusterResult[[conditionName]]$foldsClusterAssignments[[nAllFolds]],data=citrus.preclusterResult[[conditionName]]$citrus.dataArray$data,clusterCols=citrus.preclusterResult[[conditionName]]$clusterColumns))
      
      ### CHECK TO SEE IF THIS MAKES MARKERPLOT NAMES CORRECTLY 
      colLabels = citrus.preclusterResult[[conditionName]]$citrus.dataArray$fileChannelNames[[conditionName]][[1]]
      reagentNames = citrus.preclusterResult[[conditionName]]$citrus.dataArray$fileReagentNames[[conditionName]][[1]]
      displayNames = colLabels
      displayNames[nchar(reagentNames)>1] = reagentNames[nchar(reagentNames)>1]
      if (is.numeric(citrus.preclusterResult[[conditionName]]$clusterColumns)){
        colnames(clusterMedians)=displayNames[citrus.preclusterResult[[conditionName]]$clusterColumns]  
      } else {
        colnames(clusterMedians)=displayNames[colLabels%in%citrus.preclusterResult[[conditionName]]$clusterColumns]
      }
      ###
      
      rownames(clusterMedians) = conditionFeatureList[[conditionName]]$foldLargeEnoughClusters[[nAllFolds]]
      for (modelType in names(citrus.regressionResult[[conditionName]]$differentialFeatures)){
        citrus.plotHierarchicalClusterMedians(outputFile=file.path(conditionOutputDir,"markerPlots.pdf"),clusterMedians,graph=g,layout=l,plotSize=plotSize,theme=theme)
        citrus.plotHierarchicalClusterMedians(outputFile=file.path(conditionOutputDir,"markerPlotsAll.pdf"),clusterMedians,graph=g,layout=l,plotSize=plotSize,theme=theme,singlePDF=T)
        write.csv(clusterMedians,file=file.path(conditionOutputDir,"clusterMarkerMedianValues.csv"),quote=F)
        for (cvPoint in names(citrus.regressionResult[[conditionName]]$differentialFeatures[[modelType]])){
          featureClusterMatrix = .getClusterFeatureMatrix(citrus.regressionResult[[conditionName]]$differentialFeatures[[modelType]][[cvPoint]][["features"]])
          citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(conditionOutputDir,paste(modelType,"results",sep="_"),paste("featurePlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,graph=g,layout=l,plotSize=plotSize,theme=theme)
        }
      }
    }
  }
}

.graphColorPalette=function(x,alpha=1){
  #rainbow(x,alpha=.8,start=.65,end=.15)
  topo.colors(x,alpha=alpha)
}
