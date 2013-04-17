rm(list = ls())
library("citrus")

Rclusterpp.setThreads(1)

# Example 1: Diseased patients have a differing proportion of cells in clusters 2 & 3.
dataDir = file.path(system.file(package="citrus"),"extdata","example1")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c(1:2)
fileSampleSize=100

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T))
#debug(citrus.full)
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,nFolds=5,family="classification",featureTypes=c("densities"))

# Eample 2: Diseased patients have high levels of functional marker 2 in cluster 3 when stimulated. 
# No results should be visible from the unstim files.
rm(list=ls(all=T))
dataDir = file.path(system.file(package="citrus"),"extdata","example3")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")
labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"))
fileSampleSize=100
featureTypes=c("densities","medians","emDists")

conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
preclusterResult=citrus.preCluster(dataDir,outputDir,clusterCols,fileSampleSize,fileList[,-3],nFolds=5,conditionComparaMatrix=conditionComparaMatrix)
endpointResult=citrus.endpointRegress(preclusterResult,outputDir,family="classification",labels=labels,plot=T,returnResults=T,featureTypes=featureTypes,medianColumns=medianCols,emdColumns=medianCols)


# Eample 3: Diseased patients have high levels of functional marker 2 in cluster 3 when stimulated. 
# No results should be visible from the unstim files. Suvival case
rm(list=ls(all=T))
dataDir = file.path(system.file(package="citrus"),"extdata","example3")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")
labels = data.frame(time=c(floor(runif(n=10,min=11,max=20)),floor(runif(n=10,min=1,max=10))),event=sample(c(0,1,1),size=20,replace=T))
labels = data.frame(time=c(floor(runif(n=10,min=11,max=20)),floor(runif(n=10,min=1,max=10))),event=1)
nFolds=10
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"))
fileSampleSize=500
featureTypes=c("emDists")

conditionComparaMatrix=matrix(F,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[3]=T
preclusterResult=citrus.preCluster(dataDir,outputDir,clusterCols,fileSampleSize,fileList[,-3],nFolds=nFolds,conditionComparaMatrix=conditionComparaMatrix)
endpointResult=citrus.endpointRegress(preclusterResult,outputDir,family="survival",labels=labels,plot=T,returnResults=T,featureTypes=featureTypes,medianColumns=medianCols,emdColumns=medianCols)

