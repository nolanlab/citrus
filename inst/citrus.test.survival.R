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
minimumClusterSizePercent=0.1
fileList = data.frame(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),time=as.numeric(c(ceiling(c(runif(10,min=10,max=20),runif(10,min=0,max=10))))),event=as.numeric(sample(c(0,1,1),20,replace=T)))
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
featureTypes=c("densities","medians")
transformCols=NULL
family="survival"
modelTypes="glmnet"
outputDir="~/Desktop/notime/citrusTestRun/"
#citrus.quick(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
