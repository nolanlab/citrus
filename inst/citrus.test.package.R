rm(list = ls())
library("citrus")

Rclusterpp.setThreads(1)
dataDir = "Desktop/notime/citrusTestRun/"
outputDir = "Desktop/notime/citrusTestRun/citrusOutput/"
clusterCols = c(1:2)
conditionSampleSize=c(1000)
transformCols = c(1,2)

labels = rep("healthy",20)
labels[seq(1,20,by=2)]="diseased"
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T),labels=labels)

#debug(citrus.full)
citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=500,fileList=fileList,nFolds=5,featureTypes=c("medians"),medianColumns=c("Red","Blue"))
