citrus.plotTypeErrorRate = function(modelType,outputDir,regularizationThresholds,thresholdErrorRates,foldModels,cvMinima,thresholdSEMs,thresholdFDRRates=NULL){  
    nAllFolds = length(foldModels[[modelType]])
    print(paste("Plotting results for model type",modelType))
    modelOutputDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    dir.create(modelOutputDir)
    pdf(file.path(modelOutputDir,"ModelErrorRate.pdf"),width=6,height=6)
    thresholds=regularizationThresholds[[modelType]]
    errorRates=thresholdErrorRates[[modelType]]
    if (modelType=="glmnet"){
      thresholds = log(thresholds)
      xlab="log(Regularization Threshold)"
      nonzeroCounts = foldModels[[modelType]][[nAllFolds]]$df
    } 
    if (modelType=="pamr"){
      xlab="Regularization Threshold"
      nonzeroCounts = foldModels[[modelType]][[nAllFolds]]$nonzero
    }
    plot(errorRates,type='o',pch=20,col="red",main="Number of model features\n",axes=F,xlab=xlab,ylim=c(0,1),ylab="Model Cross Validation Error Rate")
    #Plot SEM
    for (i in 1:length(thresholds)){
      lines(c(i,i),c(errorRates[i]+thresholdSEMs[[modelType]][i],errorRates[i]-thresholdSEMs[[modelType]][i]),col="red",lty=3)
    }
    grid()
    axis(1,at=1:length(errorRates),labels=sapply(thresholds,citrus.formatDecimal))
    axis(2,at=c(0,.25,.5,.75,1),labels=c(0,25,50,75,100))
    axis(3,labels=nonzeroCounts,at=1:length(errorRates))
    legendLabels = c("Cross Validation Error Rate")
    legendColors = c("red")
    legendPchs = c(20)
    legendLty = c(1)
    legendPtCex=c(1)
    if (!is.null(thresholdFDRRates)&&(modelType %in% names(thresholdFDRRates))){
      lines(thresholdFDRRates[[modelType]],type='o',pch=1,cex=1.5,col="blue")
      legendLabels = c(legendLabels,"Feature False Discovery Rate")
      legendColors = c(legendColors,"blue")
      legendPchs = c(legendPchs,1)
      legendLty = c(legendLty,1)
      legendPtCex=c(legendPtCex,1)
    }
    
    cv.min = cvMinima[[modelType]]$cv.min
    cv.1se = cvMinima[[modelType]]$cv.1se
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
      
      cv.fdr.constrained = cvMinima[[modelType]]$cv.fdr.constrained
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="yellow",pch=17,cex=1.5)
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="black",pch=2,cex=1.5)
      legendLabels = c(legendLabels,"cv.fdr.constrained")
      legendColors = c(legendColors,"yellow")
      legendPchs = c(legendPchs,17)
      legendLty = c(legendLty,0)
      legendPtCex=c(legendPtCex,1.5)  
      
    }
    
    legend(x=1,y=1,legendLabels,col=legendColors,pch=legendPchs,lty=legendLty,pt.cex=legendPtCex,cex=.8,bg="white")
    
    dev.off()   
}

citrus.plotDifferentialFeatures = function(modelType,differentialFeatures,foldFeatures,outputDir,labels){
  features = foldFeatures[[length(foldFeatures)]]
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
  p = ggplot(combined) + geom_density(aes(x=value,fill=dplot,colour=dplot),alpha=.6) + facet_grid(clusterId~marker)
  print(p)
}

citrus.plotClusters = function(modelType,differentialFeatures,outputDir,clusterChildren,citrus.dataArray,conditions,clusterCols){
  data = citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions],]
  clusterChildren = clusterChildren[[length(clusterChildren)]]
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
    }
    if (nrow(data)>5000){
      bgData = data[sample(1:nrow(data),5000),clusterCols]
    } else {
      bgData = data[,clusterCols]
    }
    citrus.overlapDensityPlot(clusterDataList=clusterDataList,backgroundData=bgData)
    dev.off()
  }
}
