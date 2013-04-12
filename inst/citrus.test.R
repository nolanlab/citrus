rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
library("spade",lib.loc="~/Desktop/work/R/spade_with_nn/lib")
library("citrus")

Rclusterpp.setThreads(1)

dataDir = "inst/extdata/example3/"
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianColumns = c("FunctionalMarker1","FunctionalMarker2")
fileSampleSize=500
#transformCols = c(1,2)
nFolds=5
sapply(list.files("Desktop/work/citrus/R/",pattern=".R",full.names=T),source)

fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),class=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
featureTypes=c("densities","medians")
minimumClusterSizePercent=0.1
transformCols=NULL
family="classification"
modelTypes=c("glmnet","pamr")

citrus.quick(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
