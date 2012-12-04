rm(list = ls())
library("citrus",lib.loc="work/tmp2/")

dataDir = "work/citrus/data/syntheticData/train/unstim/"
outputDir = "work/citrus/data/testOutput/"
clusterCols = c(1:2)
conditionSampleSize=c(1000)
transformCols = c(1,2)

labels = rep("healthy",20)
labels[seq(1,20,by=2)]="diseased"
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T),labels=labels)

#debug(citrus.full)
citrus.full(dataDir,outputDir,clusterCols,conditionSampleSize,transformCols,fileList,nFolds)
