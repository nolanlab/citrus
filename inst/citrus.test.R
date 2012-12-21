rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
library("spade",lib.loc="work/R/spade_with_nn/lib")
library("citrus")

Rclusterpp.setThreads(1)

dataDir = "Desktop/work/citrus/data/syntheticData/example2"
outputDir = file.path(dataDir,"output")
clusterCols = c(1:2)
fileSampleSize=500
transformCols = c(1,2)
nFolds=5
sapply(list.files("Desktop/work/citrus/R/",pattern=".R",full.names=T),source)

fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),labels=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F


citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds)

