rm(list = ls())
library("citrus")

# Example 1: Diseased patients have a differing proportion of cells in clusters 2 & 3.
#dataDir = file.path(system.file(package="citrus"),"extdata","example1")
dataDir = file.path("~/Desktop/work/citrus/inst/extdata/example1/")
outputDir = "~/Desktop/notime/citrusOutput/"
clusterCols = c(1:2)
fileSampleSize=1000

labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
nFolds="all"
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T))
featureTypes=c("abundances")
modelTypes=c("pamr","glmnet","sam")
#modelTypes=("sam")
minimumClusterSizePercent=0.05
plot=T
returnResults=T
family="classification"
transformCols=NULL
conditionComparaMatrix=NULL
transformFactor=NULL

#res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,fileList,labels=labels,nFolds=5,family="classification",featureTypes=c("abundances"),modelTypes=c("pamr","glmnet"))
#outfile="~/Desktop/notime/glmnetprofile.out"
#Rprof(outfile)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,labels=labels,nFolds="all",family="classification",fileList=fileList,modelTypes=modelTypes,featureTypes=featureTypes,minimumClusterSizePercent=minimumClusterSizePercent,fileSampleSize=fileSampleSize,plot=T,mc.cores=4)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,labels=labels,nFolds=4,family="classification",fileList=fileList,modelTypes=modelTypes,featureTypes=featureTypes,minimumClusterSizePercent=minimumClusterSizePercent,fileSampleSize=fileSampleSize,plot=T,mc.cores=4)
#Rprof(NULL)
#summaryRprof(outfile)
citrus.endpointRegress(conditionFeatureList=res$featureObject,family="classification",modelTypes="pamr",labels=labels)
filePopulationList = list(unstim=matrix(list.files("~/Desktop/work/citrus/inst/extdata/example4.1",pattern=".fcs"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))))
res = citrus.full(dataDir,outputDir,clusterCols,fileSampleSize,filePopulationList=filePopulationList,labels=labels,family="classification",featureTypes=c("abundances"),modelTypes=c("glmnet","pamr"),plot=T)


# Example 2: Diseased patients have high levels of functional marker 2 in cluster 3 when stimulated. 
# No results should be visible from the unstim files.
rm(list=ls(all=T))
dataDir = file.path(system.file(package="citrus"),"extdata","example3")
#dataDir = file.path("~/Desktop/work/citrus/inst/extdata/example3/")
outputDir = "~/Desktop/notime/tmp/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")
labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"))
fileSampleSize=1000
featureTypes=c("abundances","medians","emDists")
featureTypes=c("emDists")
family="classification"
modelTypes="glmnet"
modelTypes=c("pamr","glmnet")
minimumClusterSizePercent=0.05
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F

# Results should be bad
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds=5,family=family,fileList=fileList,modelTypes=modelTypes,featureTypes=c("abundances"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix)
# Results should be pretty good
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds=5,family=family,fileList=fileList,modelTypes=modelTypes,featureTypes=c("medians"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,medianColumns=medianCols)
# Results should be pretty good
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds=5,family=family,fileList=fileList,modelTypes=modelTypes,featureTypes=c("emDists"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,emdColumns=medianCols)

# Same but quicker
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds="all",family="classification",fileList=fileList,modelTypes=c("glmnet","pamr"),featureTypes=c("abundances"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds="all",family="classification",fileList=fileList,modelTypes=c("glmnet","pamr"),featureTypes=c("medians"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,medianColumns=medianCols)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,nFolds="all",family="classification",fileList=fileList,modelTypes=c("glmnet","pamr"),featureTypes=c("emDists"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,emdColumns=medianCols)

# Same but using hand-gated data
dataDir = file.path(system.file(package="citrus"),"extdata","example2.1")
filePopulationList = list(unstim=matrix(list.files(dataDir,pattern="unstim"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))),stim1=matrix(list.files(dataDir,pattern="stim1"),ncol=3,byrow=T,dimnames=list(NULL,paste("Pop",1:3))))
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,family="classification",filePopulationList=filePopulationList,modelTypes=c("glmnet","pamr"),featureTypes=c("abundances"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,family="classification",filePopulationList=filePopulationList,modelTypes=c("glmnet","pamr"),featureTypes=c("medians"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,medianColumns=medianCols)
res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=fileSampleSize,labels=labels,family="classification",filePopulationList=filePopulationList,modelTypes=c("glmnet","pamr"),featureTypes=c("emDists"),minimumClusterSizePercent=0.05,conditionComparaMatrix=conditionComparaMatrix,emdColumns=medianCols)

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

# Example 4: multiclass
# Class 1 has small amounts of population 3 (1%), class 2 has small more of population 2 (%), and class 3 has lots of population 2 and 3 (8%)
rm(list = ls())
library("citrus")
dataDir = file.path("~/Desktop/work/citrus/inst/extdata/example4/")
outputDir = "~/Desktop/notime/citrusOutput/"
clusterCols = c(1:2)
fileSampleSize=500

labels = c(rep("Class1",10),rep("Class2",10),rep("Class3",10))
labels = as.factor(labels)
nFolds=5
nFolds="all"
fileList = data.frame(unstim=list.files(dataDir,pattern=".fcs",ignore.case=T))
featureTypes=c("abundances")
#modelTypes=c("pamr","glmnet","sam")
modelTypes=("glmnet")
minimumClusterSizePercent=0.05
plot=T
returnResults=T
family="classification"
transformCols=NULL
conditionComparaMatrix=NULL
transformFactor=NULL

res = citrus.full(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,labels=labels,nFolds="all",family=family,fileList=fileList,modelTypes=modelTypes,featureTypes=featureTypes,minimumClusterSizePercent=minimumClusterSizePercent,fileSampleSize=fileSampleSize,plot=T,mc.cores=4)




# Testing the new file mapping example
rm(list=ls(all=T))
dataDir = file.path(system.file(package="citrus"),"extdata","example3")
#dataDir = file.path("~/Desktop/work/citrus/inst/extdata/example3/")
outputDir = "~/Desktop/notime/citrusOutput/"
clusterCols = c("LineageMarker1","LineageMarker2")
medianCols = c("FunctionalMarker1","FunctionalMarker2")
labels = c(rep("healthy",10),rep("diseased",10))
labels = as.factor(labels)
nFolds=5
fileList = cbind(unstim=list.files(dataDir,pattern="unstim"),stim1=list.files(dataDir,pattern="stim1"))
fileSampleSize=1000
featureTypes=c("abundances","medians","emDists")
featureTypes=c("emDists")
family="classification"
modelTypes="glmnet"
modelTypes="sam"
minimumClusterSizePercent=0.05
conditionComparaMatrix=matrix(T,ncol=2,nrow=2,dimnames=list(c("unstim","stim1"),c("unstim","stim1")))
conditionComparaMatrix[2]=F

trainFileList = fileList[seq(from=1,to=19,by=2),]
testFileList = fileList[seq(from=2,to=20,by=2),]
preClusterResult = citrus.preCluster(dataDir=dataDir,outputDir=outputDir,clusterCols=clusterCols,fileSampleSize=1000,fileList=trainFileList,nFolds="all",conditionComparaMatrix=conditionComparaMatrix,mc.cores=4)
mappingResults = citrus.mapFileDataToClustering(dataDir=dataDir,newFileList=testFileList,fileSampleSize=1000,preClusterResult=preClusterResult,mc.cores=4)

trainFeatures = citrus.buildFeatures(preclusterResult=preClusterResult,outputDir=outputDir,featureTypes=c("medians"),minimumClusterSizePercent=0.01,medianColumns=medianCols)
trainLargeEnoughClusters = lapply(names(trainFeatures),.extractConditionLargeEnoughClusters,foldFeatures=trainFeatures)
names(trainLargeEnoughClusters) = names(trainFeatures)
trainFeatures

citrus.exportConditionClusters(conditionClusterIds=list(unstim_vs_stim1=c(15800,12222),stim1=c(8839,9987)),preclusterResult=preClusterResult,outputDir=outputDir,sampleIds=c(6,10))
f = trainFeatures$unstim_vs_stim1$foldFeatures[[1]]
cvm = cv.glmnet(x=f,y=as.factor(rep(c("H","D"),each=5)),family="binomial",type.measure="class")
plot(cvm)
fastICA(X=f,n.comp=2)
cvm = cv.glmnet(x=fastICA(f,n.comp=3)$S,y=as.factor(rep(c("H","D"),each=5)),family="binomial",type.measure="class")
plot(cvm)




