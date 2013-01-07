rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
library("spade",lib.loc="~/Desktop/work/R/spade_with_nn/lib")
library("citrus")

Rclusterpp.setThreads(1)

dataDir = "Desktop/work/citrus/data/syntheticData/example3"
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianColumns = c("FunctionalMarker1","FunctionalMarker2")
fileSampleSize=1000
#transformCols = c(1,2)
nFolds=5
sapply(list.files("Desktop/work/citrus/R/",pattern=".R",full.names=T),source)

fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),labels=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(F,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[3]=T
featureTypes=c("densities","medians")


citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))

