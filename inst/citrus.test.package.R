rm(list = ls())
library("citrus")

Rclusterpp.setThreads(1)

# Example 1: Diseased patients have a differing proportion of cells in clusters 2 & 3.
dataDir = file.path(system.file(package="citrus"),"data","example1")
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c(1:2)
conditionSampleSize=c(1000)

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T),labels=labels)

#debug(citrus.full)
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=500,fileList=fileList,nFolds=5,featureTypes=c("densities"))


# Eample 2: Diseased patients have high levels of functional marker 2 in cluster 3 when stimulated. 
# No results should be visible from the unstim files.
rm(list=ls(all=T))
dataDir = "Desktop/work/citrus/data/syntheticData/example2/"
outputDir = "Desktop/notime/citrusTestRun/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),labels=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
fileSampleSize=1000

# First just calculate density features. There should not be any good results here.
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,featureTypes=c("densities"))
# Now try with medians. Should get good results.
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,featureTypes=c("medians"),medianColumns=c("FunctionalMarker1","FunctionalMarker2"))
# Now try with both medians and densities. Should get good results and automatically find the good features.
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,featureTypes=c("densities","medians"),medianColumns=c("FunctionalMarker1","FunctionalMarker2"))

fileSampleSize=1000
#transformCols = c(1,2)
nFolds=5
sapply(list.files("Desktop/work/citrus/R/",pattern=".R",full.names=T),source)

fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"),labels=c(rep("healthy",10),rep("diseased",10)))
conditionComparaMatrix=matrix(F,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[3]=T
featureTypes=c("densities","medians")


citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,fileList=fileList,nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix,featureTypes="medians",medianColumns=c("FunctionalMarker1","FunctionalMarker2"))

