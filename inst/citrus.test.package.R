rm(list = ls())
library("citrus")

#Rclusterpp.setThreads(1)

# Example 1: Diseased patients have a differing proportion of cells in clusters 2 & 3.
dataDir = file.path(system.file(package="citrus"),"extdata","example1")
outputDir = "~/Desktop/notime/tmp/citrusOutput/"
clusterCols = c(1:2)
fileSampleSize=1000

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T))
res = citrus.readFCSSet(dataDir=dataDir,fileList=fileList,conditions="unstim",fileSampleSize=50)
featureTypes=c("densities")
modelTypes=c("pamr","glmnet")
minimumClusterSizePercent=0.05
plot=T
returnResults=T
family="classification"
transformCols=NULL
conditionComparaMatrix=NULL
transformFactor=NULL
#Sys.setenv(OMP_NUM_THREADS=1)
nFolds=5
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels=labels,nFolds=5,family="classification",featureTypes=c("densities"))
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels=labels,nFolds="all",family="classification",featureTypes=c("densities"),plot=T)


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

f=endpointResult$unstim_vs_stim1$foldFeatures[[11]]
dim(f)
labels = data.frame(time=c(floor(runif(n=10,min=20,max=70)),floor(runif(n=10,min=1,max=10))),event=1)
s=Surv(time=labels[,1],event=labels[,2])
plot(f[,54])
plot(f[,54],labels[,1])
s=Surv(time=((f[,42]+2)*5)+rnorm(20),event=labels[,2])
plot(f[,42],s[,1])
labels = data.frame(time=s[,1],event=s[,2])
plot(glmnet(x=res$unstim_vs_stim1$features,y=s,family="cox"))
plot(cv.glmnet(x=f,y=s,family="cox"))
sm=SAM(x=t(f),y=s[,1],censoring.status=s[,2],resp.type="Survival",fdr.output=0.000,genenames=colnames(f),nperms=10000)
sm
df = data.frame(f)
coxph(s~cluster.19981.FunctionalMarker2.emDist,data=df)
