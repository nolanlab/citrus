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
wd = "work/citrus/testOutput/"

clusterCols = c(1:2)
fileNames = cbind(unstim=list.files(dataDir,pattern="train"))

citrus.dataArray = citrus.readFCSSet(dataDirectory=dataDir,conditionFileList=fileNames,conditionSampleSize=c(500),transformColumns=c(1,2))

labels = rep("healthy",20)
labels[grep("diseased",fileNames[,1])]="diseased"
labels = as.factor(labels)

folds = citrus.leavoutFold(x=1,y=labels,leaveoutSize=5)
folds = pamr:::balanced.folds(y=labels,nfolds=5)
folds[[(length(folds)+1)]]="all"

#clusterConditionMatrix = matrix(c(T,F,T,T),ncol=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
clusterConditionMatrix = matrix(c(T),ncol=1,dimnames=list(c("unstim"),c("unstim")))

conditions="unstim"
foldsCluster = lapply(folds,citrus.foldCluster,citrus.dataArray=citrus.dataArray,clusterCols=c(1,2),conditions="unstim")
foldsClusterAssignments = lapply(foldsCluster,citrus.calculateCompleteHierarchicalMembership)
leftoutClusterAssignments = lapply(1:5,citrus.mapFoldDataToClusterSpace,citrus.dataArray=citrus.dataArray,foldClusterAssignments=foldsClusterAssignments,folds=folds,conditions=conditions,clusterColumns=c(1,2))
foldLargeEnoughClusters = lapply(1:6,citrus.calculateFoldLargeEnoughClusters,foldsClusterAssignments=foldsClusterAssignments,folds=folds,citrus.dataArray=citrus.dataArray)

#conditions="unstim"
foldFeatures = lapply(1:6,citrus.buildFoldFeatures,featureTypes=c("densities","medians"),citrus.dataArray=citrus.dataArray,foldsClusterAssignments=foldsClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=c(1,2))
leftoutFeatures = lapply(1:5,citrus.buildFoldFeatures,featureTypes=c("densities","medians"),citrus.dataArray=citrus.dataArray,foldsClusterAssignments=leftoutClusterAssignments,foldLargeEnoughClusters=foldLargeEnoughClusters,conditions=conditions,medianColumns=c(1,2),calculateLeaveoutData=T)

modelTypes = c("pamr","glmnet")
regularizationThresholds = citrus.generateRgularizationThresholds(foldFeatures[[6]],labels,modelTypes=modelTypes)

thresholdErrorRates = list()
finalModels = list();
for (modelType in modelTypes){
  foldModels = lapply(1:5,citrus.buildFoldModels,folds,foldFeatures,labels,type=modelType,regularizationThresholds[[modelType]])
  finalModels[[modelType]] = citrus.buildModel(features=foldFeatures[[6]],labels=labels,type=modelType,regularizationThresholds=regularizationThresholds[[modelType]])
  leftoutPredictions = lapply(1:5,citrus.foldPredict,model=foldModels,features=leftoutFeatures,regularizationThresholds=regularizationThresholds)
  predictionSuccess = lapply(1:5,citrus.foldScore,folds=folds,predictions=leftoutPredictions,labels=labels)
  thresholdErrorRates[[modelType]] = 1-(apply(do.call("rbind",predictionSuccess),2,sum)/nrow(do.call("rbind",predictionSuccess)))
}

pamrFdrCount = pamr.fdr(finalModels[["pamr"]],data=list(x=t(foldFeatures[[6]]),y=labels),nperms=1000)

plot(regularizationThresholds$glmnet,thresholdErrorRates[["glmnet"]],ylim=c(0,.5),type='o')
plot(regularizationThresholds$pamr,thresholdErrorRates[["pamr"]],ylim=c(0,.5),type='o',col="red")
lines(pamrFdrCount$results[,"Threshold"],pamrFdrCount$results[,"Median FDR"],type='o',col="blue")
