citrus.plotTypeErrorRate = function(modelType,outputDir,regularizationThresholds,thresholdErrorRates,foldModels,cvMinima,thresholdSEMs,thresholdFDRRates=NULL){  
    nAllFolds = length(foldModels[[modelType]])
    print(paste("Plotting results for model type",modelType))
    modelOutputDir = file.path(outputDir,paste(modelType,"_results/",sep=""))
    dir.create(modelOutputDir)
    pdf(file.path(modelOutputDir,"ModelErrorRate.pdf"),width=8,height=8)
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
    plot(errorRates,type='o',pch=20,col="red",main=paste(modelType,"model error rate\n"),axes=F,xlab=xlab,ylim=c(0,1),ylab="Percent")
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
    if (modelType %in% names(thresholdFDRRates)){
      lines(thresholdFDRRates[[modelType]],type='o',pch=1,cex=1.5,col="blue")
      legendLabels = c(legendLabels,"Feature False Discovery Rate")
      legendColors = c(legendColors,"blue")
      legendPchs = c(legendPchs,1)
      legendLty = c(legendLty,1)
      legendPtCex=c(legendPtCex,1)
    }
    
    cv.min = cvMinima[[modelType]]$cv.min
    cv.1se = cvMinima[[modelType]]$cv.1se
    points(c(cv.min,cv.min),y=c(errorRates[cv.min],errorRates[cv.min]),col="green",pch=20,cex=2)
    points(c(cv.1se,cv.1se),y=c(errorRates[cv.1se],errorRates[cv.1se]),col="orange",pch=9,cex=2)
    
    
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
    modelTypeDir = paste(outputDir,modelType,"_results/",sep="")
    nonzeroFeatureNames = differentialFeatures[[modelType]][[cvPoint]][["features"]]
    for (nonzeroFeatureName in nonzeroFeatureNames){
      df = data.frame(value=features[,nonzeroFeatureName],labels,featureName=nonzeroFeatureName)
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
    p = p + facet_wrap(~featureName,ncol=4) + geom_boxplot(outlier.colour=rgb(0,0,0,0),colour=rgb(0,0,0,.3)) + geom_point(aes(color=labels),alpha=I(0.25),shape=19,size=I(2)) + scale_colour_manual(values = c("red","blue")) + coord_flip() +  theme_bw() + ylab("") + xlab("") + opts(legend.position = "none") 
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

citrus.overlapDensityPlot = function(data,backgroundRef=NULL,scaleUp=.8){
  samples = names(data)
  ncol=ncol(data[[samples[1]]])
  
  nrow=length(names(data))
  
  dMax = max(unlist(lapply(data,max)))
  dMin = min(unlist(lapply(data,min)))
  if (!is.null(backgroundRef)){
    dAll = apply(backgroundRef,2,density,from=dMin,to=dMax)
  }
  par(mar=c(1,5,0,0),oma=c(0,0,0,0))
  
  
  plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(1,ncol+1),ylim=c(1,nrow+1+scaleUp))
  axis(side=2,labels=names(data),at=c(1:nrow),lwd=0,lwd.ticks=1,las=1,line=-1,cex.axis=1.3)
  #text(x=ncol/2,y=0,labels=paste("Range: ",sprintf("%1.1f",dMin),"-",sprintf("%1.1f",dMax)))
  for (i in 1:nrow){
    for (j in 1:ncol){
      d = density(data[[samples[i]]][,j],from=dMin,to=dMax)
      range = scaleToRange(dAll[[j]]$x,scale=c(j,j+.8))
      x0 = range[which(abs(dAll[[j]]$x)==(min(abs(dAll[[j]]$x))))]
      xMin = range[which(dAll[[j]]$x==min(dAll[[j]]$x))]
      xMax = range[which(dAll[[j]]$x==max(dAll[[j]]$x))]
      
      lines(x=c(x0,x0),y=c(0,(nrow+scaleUp)),col=rgb(0,0,1,0.3),lty=3,lwd=1)
      
      if (!is.null(backgroundRef)){
        lines(scaleToRange(dAll[[j]]$x,scale=c(j,(j+.8))),scaleToRange(dAll[[j]]$y,scale=c(i,(i+scaleUp))),type='l',col=rgb(0,0,0,.3),lty=1,lwd=2)
        lines(scaleToRange(dAll[[j]]$x,scale=c(j,(j+.8))),scaleToRange(dAll[[j]]$y,scale=c(i,(i+scaleUp))),type='l',col=rgb(0,0,0,1),lty=5,lwd=2)
        
      }
      lines(scaleToRange(d$x,scale=c(j,(j+.8))),scaleToRange(d$y,scale=c(i,(i+scaleUp))),type='l',col="red",lwd=2)
      
    }
  }
  for (j in 1:ncol){
    range = scaleToRange(dAll[[j]]$x,scale=c(j,j+.8))
    x0 = range[which(abs(dAll[[j]]$x)==(min(abs(dAll[[j]]$x))))]
    xMin = range[which(dAll[[j]]$x==min(dAll[[j]]$x))]
    xMax = range[which(dAll[[j]]$x==max(dAll[[j]]$x))]
    
    #lines(x=c(x0,x0),y=c(0,(nrow+scaleUp)),col=rgb(0,0,1,0.5),lty=1,lwd=2)
    lines(x=c(xMin,xMin),y=c(0,(nrow+scaleUp)),col="grey",lty=2,lwd=2)
    lines(x=c(xMax,xMax),y=c(0,(nrow+scaleUp)),col="grey",lty=2,lwd=2)
  }
  
  text(x=(1:ncol)-.2,y=nrow+scaleUp+.4,labels=colnames(data[[samples[1]]]),cex=1.2,pos=4)
  #text(x=(1:ncol)+.2,y=nrow+scaleUp+.4,labels=colnames(data[[samples[1]]]))
  print(paste("Range: ",sprintf("%1.1f",dMin),"-",sprintf("%1.1f",dMax)))
}


citrus.plotClusters = function(modelType,differentialFeatures,outputDir,clusterChildren,citrus.dataArray,conditions,clusterCols){
  data = citrus.dataArray$data[citrus.dataArray$data[,"fileId"] %in% citrus.dataArray$fileIds[,conditions],]
  clusterChildren = clusterChildren[[length(clusterChildren)]]
  for (cvPoint in names(differentialFeatures[[modelType]])){
    nonzeroClusters = as.numeric(differentialFeatures[[modelType]][[cvPoint]][["clusters"]])
    pdf(file=file.path(outputDir,paste(modelType,"_results/clusters-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=(2*length(clusterCols)+2),height=(2*length(nonzeroClusters)+2))
    clusterDataList=list();
    for (nonzeroCluster in sort(nonzeroClusters)){
      clusterDataList[[as.character(nonzeroCluster)]]=data[clusterChildren[[nonzeroCluster]],clusterCols]
    }
    citrus.overlapDensityPlot(clusterDataList,backgroundRef=data[,clusterCols])
    dev.off()
  }
}