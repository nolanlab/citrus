rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("spade",lib.loc="work/R/spade_with_nn/lib")

source("work/citrus/citrus.cluster.R")
source("work/citrus/citrus.util.R")
source("work/citrus/citrus.featureFunctions.R")
source("work/citrus/citrus.classificationModel.R")

dataDir = "work/citrus/syntheticData/unstim/"
outputDir = "work/citrus/testOutput/"

clusterCols = c(1:2)
fileNames = cbind(unstim=list.files(dataDir,pattern="train"))

citrus.dataArray = citrus.readFCSSet(dataDirectory=dataDir,conditionFileList=fileNames,conditionSampleSize=c(500),transformColumns=c(1,2))

labels = rep("healthy",20)
labels[grep("diseased",fileNames[,1])]="diseased"
labels = as.factor(labels)

nFolds=5
nAllFolds = nFolds+1
folds = citrus.leavoutFold(x=1,y=labels,leaveoutSize=5)
folds = pamr:::balanced.folds(y=labels,nfolds=nAllFolds)
folds[[nAllFolds]]="all"

#clusterConditionMatrix = matrix(c(T,F,T,T),ncol=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
clusterConditionMatrix = matrix(c(T),ncol=1,dimnames=list(c("unstim"),c("unstim")))

conditions="unstim"
foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=c(1,2),conditions="unstim")
foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)
leftoutClusterAssignments = lapply(1:nFolds,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterColumns=c(1,2))
foldLargeEnoughClusters = lapply(1:nAllFolds,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=foldsClusterAssignments,folds=folds,citrus.dataArray=citrus.dataArray)

foldFeatures = lapply(1:nAllFolds,citrus.buildFoldFeatures,featureTypes=c("densities"),citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=c(1,2))
leftoutFeatures = lapply(1:nFolds,citrus.buildFoldFeatures,featureTypes=c("densities"),citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=c(1,2),calculateLeaveoutData=T)

modelTypes = c("glmnet","pamr")
regularizationThresholds = citrus.generateRgularizationThresholds(foldFeatures[[nAllFolds]],labels,modelTypes=modelTypes)

log(regularizationThresholds$glmnet)
plot(cv.glmnet(x=foldFeatures[[nAllFolds]],y=labels,family="binomial",type.measure="class"))
pd = list(x=t(foldFeatures[[nAllFolds]]),y=labels)
pm = pamr.train(data=pd)
cv = pamr.cv(pm,pd)
pamr.plotcv(cv)

thresholdErrorRates = list();
thresholdFDRRates = list();
thresholdSEMs = list();
foldModels = list();
for (modelType in modelTypes){
  foldModels[[modelType]] = lapply(1:nAllFolds,citrus.buildFoldModels,folds=folds,foldFeatures=foldFeatures,labels=labels,type=modelType,regularizationThresholds=regularizationThresholds[[modelType]])
  leftoutPredictions = lapply(1:nFolds,citrus.foldPredict,model=foldModels[[modelType]],features=leftoutFeatures)
  predictionSuccess = lapply(1:nFolds,citrus.foldScore,folds=folds,predictions=leftoutPredictions,labels=labels)
  thresholdSEMs[[modelType]] = citrus.calculateSEM(predictionSuccess)
  thresholdErrorRates[[modelType]] = 1-(apply(do.call("rbind",predictionSuccess),2,sum)/nrow(do.call("rbind",predictionSuccess)))
  if (modelType=="pamr"){
    thresholdFDRRates[[modelType]] = pamr.fdr.new(foldModels[[modelType]][[nAllFolds]],data=list(x=t(foldFeatures[[nAllFolds]]),y=labels),nperms=1000)$results[,"Median FDR"]
  }    
}

for (modelType in modelTypes){
  print(paste("Plotting results for model type",modelType)  )
  modelOutputDir = paste(outputDir,modelType,"_results/",sep="")
  dir.create(modelOutputDir)
  thresholds=regularizationThresholds[[modelType]]
  errorRates=thresholdErrorRates[[modelType]]
  if (modelType=="glmnet"){
    thresholds = log(thresholds)
    xlab="log(Regularization Threshold)"
  } 
  if (modelType=="pamr"){
    xlab="Regularization Threshold"
  }
  plot(errorRates,type='o',pch=20,col="red",main=paste(modelType,"model error rate\n"),axes=F,xlab=xlab,ylim=c(0,1),ylab="Percent")
  #Plot SEM
  for (i in 1:length(thresholds)){
    lines(c(i,i),c(errorRates[i]+thresholdSEMs[[modelType]][i],errorRates[i]-thresholdSEMs[[modelType]][i]),col="red",lty=3)
  }
  grid()
  axis(1,at=1:length(errorRates),labels=sapply(thresholds,citrus.formatDecimal))
  axis(2,at=c(0,.25,.5,.75,1),labels=c(0,25,50,75,100))
  axis(3,labels=foldModels[[modelType]][[nAllFolds]]$nonzero,at=1:length(errorRates))
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
  
  cv.min = min(which(errorRates==min(errorRates)))
  cv.1se = min(which(errorRates<=(errorRates[cv.min]+thresholdSEMs[[modelType]][cv.min])))
  points(c(cv.min,cv.min),y=c(errorRates[cv.min],errorRates[cv.min]),col="green",pch=20,cex=2)
  points(c(cv.1se,cv.1se),y=c(errorRates[cv.1se],errorRates[cv.1se]),col="orange",pch=9,cex=2)
  
  
  legendLabels = c(legendLabels,"CV.min","CV.1se")
  legendColors = c(legendColors,"green","orange")
  legendPchs = c(legendPchs,20,9)
  legendLty = c(legendLty,0,0)
  legendPtCex=c(legendPtCex,2,1.5)  
  
  
  if (modelType %in% names(thresholdFDRRates)){
    cv.fdr.constrained = max(intersect(which(thresholdFDRRates[[modelType]]<0.01),which(errorRates==min(errorRates))))
    points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="yellow",pch=17,cex=1.5)
    points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="black",pch=2,cex=1.5)
    legendLabels = c(legendLabels,"CV.FDR.Constrained")
    legendColors = c(legendColors,"yellow")
    legendPchs = c(legendPchs,17)
    legendLty = c(legendLty,0)
    legendPtCex=c(legendPtCex,1.5)  
    
  }
  
  legend(x=1,y=1,legendLabels,col=legendColors,pch=legendPchs,lty=legendLty,pt.cex=legendPtCex,cex=.8)
}

