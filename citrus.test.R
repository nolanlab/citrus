rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
library("spade",lib.loc="work/R/spade_with_nn/lib")

source("work/citrus/citrus.cluster.R")
source("work/citrus/citrus.util.R")
source("work/citrus/citrus.featureFunctions.R")
source("work/citrus/citrus.classificationModel.R")
source("work/citrus/citrus.external.R")
source("work/citrus/citrus.plot.R")


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

foldModels = lapply(modelTypes,citrus.buildTypeModels,folds=folds,foldFeatures=foldFeatures,labels=labels,regularizationThresholds=regularizationThresholds)
names(foldModels)=modelTypes

leftoutPredictions = lapply(modelTypes,citrus.foldTypePredict,foldModels=foldModels,leftoutFeatures=leftoutFeatures)
names(leftoutPredictions)=modelTypes

predictionSuccess = lapply(as.list(modelTypes),citrus.foldTypeScore,folds=folds,leftoutPredictions=leftoutPredictions,labels=labels)
names(predictionSuccess)=modelTypes

thresholdSEMs = lapply(modelTypes,citrus.modelTypeSEM)
names(thresholdSEMs)=modelTypes

thresholdErrorRates = lapply(modelTypes,citrus.calcualteTypeErroRate,predictionSuccess=predictionSuccess)
names(thresholdErrorRates)=modelTypes

thresholdFDRRates = lapply(modelTypes,citrus.calculateTypeFDRRate,foldModels=foldModels,foldFeatures=foldFeatures,labels=labels)
names(thresholdFDRRates)=modelTypes

cvMinima = lapply(modelTypes,citrus.getCVMinima,thresholdErrorRates=thresholdErrorRates,thresholdSEMs=thresholdSEMs,thresholdFDRRates=thresholdFDRRates)
names(cvMinima)=modelTypes

# Plot
sapply(modelTypes,citrus.plotTypeErrorRate,outputDir=outputDir,regularizationThresholds=regularizationThresholds,thresholdErrorRates=thresholdErrorRates,thresholdFDRRates=thresholdFDRRates,cvMinima=cvMinima)

# Extract Features
differentialFeatures = lapply(modelTypes,citrus.extractModelFeatures,cvMinima=cvMinima,foldModels=foldModels,foldFeatures=foldFeatures,regularizationThresholds=regularizationThresholds)
names(differentialFeatures) = modelTypes

# Plot Features
lapply(modelTypes,citrus.plotDifferentialFeatures,differentialFeatures=differentialFeatures,foldFeatures=foldFeatures,outputDir=outputDir,labels=labels)

# Plot Clusters
lapply(modelTypes,citrus.plotClusters,differentialFeatures=differentialFeatures,outputDir=outputDir,clusterChildren=foldsClusterAssignments,citrus.dataArray=citrus.dataArray,conditions=conditions,clusterCols=clusterCols)
