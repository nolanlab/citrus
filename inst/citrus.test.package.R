rm(list = ls())
library("citrus")

# Example 1: Diseased patients have a differing proportion of cells in clusters 2 & 3.
dataDir = file.path(system.file(package="citrus"),"extdata","example1")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c(1:2)
fileSampleSize=1000

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T))
featureTypes=c("densities")
modelTypes=c("pamr","glmnet")
minimumClusterSizePercent=0.05
plot=T
returnResults=T
family="classification"
transformCols=NULL
conditionComparaMatrix=NULL
transformFactor=NULL
nFolds=5
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels=labels,nFolds=5,family="classification",featureTypes=c("densities"))
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels=labels,nFolds="all",family="classification",featureTypes=c("densities"),plot=T)

filePopulationList = list(unstim=matrix(list.files("~/Desktop/work/citrus/inst/extdata//example4.1",pattern=".fcs"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))))
filePopulationList = list(unstim=matrix(list.files("~/Desktop/work/citrus/inst/extdata//example4.1",pattern=".fcs"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))),stim1=matrix(list.files("~/Desktop/work/citrus/inst/extdata//example4.1",pattern=".fcs"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))))
x = citrus.assembleHandGates(dataDir=dataDir,filePopulationList=filePopulationList)

conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F

x = citrus.assembleHandGates(dataDir=dataDir,filePopulationList=filePopulationList,conditionComparaMatrix=conditionComparaMatrix)
length(x$unstim_vs_stim1$foldsClusterAssignments$all)
length(x$stim1$foldsClusterAssignments$all)
table(x$stim1$citrus.dataArray$data[x$stim1$foldsClusterAssignments$all[[3]],"fileId"])



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
fileSampleSize=1000
featureTypes=c("densities","medians","emDists")
family="classification"
modelTypes="glmnet"
minimumClusterSizePercent=0.05
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F
preclusterResult=citrus.preCluster(dataDir,outputDir,clusterCols,fileSampleSize,fileList[,-3],nFolds=5,conditionComparaMatrix=conditionComparaMatrix)
endpointResult=citrus.endpointRegress(preclusterResult,outputDir,family="classification",labels=labels,plot=T,returnResults=T,featureTypes=featureTypes,medianColumns=medianCols,emdColumns=medianCols)
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,nFolds=5,family,modelTypes,featureTypes=featureTypes,minimumClusterSizePercent,conditionComparaMatrix=conditionComparaMatrix,plot=T,returnResults=T,medianColumns=medianCols,emdColumns=medianCols)
res = citrus.quick(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,family,modelTypes,featureTypes=featureTypes,minimumClusterSizePercent,conditionComparaMatrix=conditionComparaMatrix,plot=T,returnResults=T,medianColumns=medianCols,emdColumns=medianCols)

conditionComparaMatrix=matrix(F,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[3]=T
res = citrus.quick(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,family,modelTypes,featureTypes="emDists",minimumClusterSizePercent,conditionComparaMatrix=conditionComparaMatrix,plot=T,returnResults=T,medianColumns=medianCols,emdColumns=medianCols)


# Eample 3: Diseased patients have high levels of functional marker 2 in cluster 3 when stimulated. 
# No results should be visible from the unstim files. Suvival case
rm(list=ls(all=T))
dataDir = file.path(system.file(package="citrus"),"extdata","example3")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")
labels = data.frame(time=c(floor(runif(n=10,min=13,max=30)),floor(runif(n=10,min=1,max=10))),event=sample(c(0,1,1),size=20,replace=T))
labels = data.frame(time=c(floor(runif(n=10,min=10,max=20)),floor(runif(n=10,min=1,max=10))),event=1)
nFolds=10
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"))
fileSampleSize=500
featureTypes=c("emDists")
family="survival"
modelTypes="glmnet"
minimumClusterSizePercent=0.1
Rclusterpp.setThreads(1)


conditionComparaMatrix=matrix(F,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[3]=T
res=citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,nFolds=5,family,modelTypes="glmnet",featureTypes="emDists",minimumClusterSizePercent=minimumClusterSizePercent,conditionComparaMatrix=conditionComparaMatrix,plot=T,returnResults=T,emdColumns=medianCols)
res=citrus.quick(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels,family,modelTypes="glmnet",featureTypes="emDists",minimumClusterSizePercent=minimumClusterSizePercent,conditionComparaMatrix=conditionComparaMatrix,plot=T,returnResults=T,emdColumns=medianCols)

