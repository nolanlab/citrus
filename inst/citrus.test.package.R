rm(list = ls())
library("citrus")

Rclusterpp.setThreads(1)
dataDir = "Desktop/work/citrus/data/syntheticData/example1/"
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c(1:2)
conditionSampleSize=c(1000)
transformCols = c(1,2)

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T),labels=labels)

#debug(citrus.full)
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=500,fileList=fileList,nFolds=5,featureTypes=c("densities"))


# Eample 3
rm(list=ls(all=T))
dataDir = "Desktop/work/citrus/data/syntheticData/example3/"
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),labels=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
fileSampleSize=500

#debug(citrus.full)
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds)
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=1000,fileList=fileList,nFolds=nFolds,featureTypes=c("medians"),medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
