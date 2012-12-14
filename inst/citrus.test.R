rm(list = ls())
library("flowCore")
library("Rclusterpp")
library("pamr")
library("glmnet")
library("ggplot2")
library("spade",lib.loc="work/R/spade_with_nn/lib")

source("Desktop/work/citrus/R/citrus.cluster.R")


dataDir = "work/citrus_old/syntheticData/train/unstim/"
outputDir = "work/citrus_old/testOutput/"
clusterCols = c(1:2)
conditionSampleSize=c(500)
transformCols = c(1,2)

labels = rep("healthy",20)
labels[seq(1,20,by=2)]="diseased"
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T),labels=labels)

citrus.full(dataDir,outputDir,clusterCols,conditionSampleSize,transformCols,fileList,nFolds)