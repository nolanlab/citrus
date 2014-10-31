########################
# Helper Functions
########################
.graphColorPalette=function(x,alpha=1){
  #rainbow(x,alpha=.8,start=.65,end=.15)
  topo.colors(x,alpha=alpha)
}

.formatDecimal = function(x){
  sprintf("%1.2f", x)
}

.getDisplayNames=function(citrus.combinedFCSSet,clusteringColumns){
  colLabels = citrus.combinedFCSSet$fileChannelNames[[1]][[1]]
  reagentNames = citrus.combinedFCSSet$fileReagentNames[[1]][[1]]
  displayNames = colLabels
  displayNames[nchar(reagentNames)>2] = reagentNames[nchar(reagentNames)>2]
  if (all(is.numeric(clusteringColumns))){
    return(displayNames[clusteringColumns])
  } else {
    names(displayNames) = colLabels
    return(as.vector(displayNames[clusteringColumns]))
  }
  
}

.getClusterMedians = function(clusterId,clusterAssignments,clusterCols,data){
  apply(data[clusterAssignments[[clusterId]],clusterCols],2,median)
}

.decimalFormat = function(x){
  sprintf("%.2f",x)
}

.scaleToOne = function(x){
  x = x-min(x)
  x = x/max(x)
  return(x)
}

.getClusterFeatureMatrix = function(featureVector){
  df = do.call("rbind",strsplit(gsub(pattern="(cluster [0-9]+) ",replacement="\\1~",featureVector),"~"))
  return(cbind(cluster=do.call("rbind",strsplit(df[,1]," "))[,2],feature=df[,2]))
}

citrus.plotTypeErrorRate = function(modelType,modelOutputDirectory,regularizationThresholds,thresholdCVRates,finalModel,cvMinima,family){  
    if (modelType=="sam"){
      return()
    }
    
    pdf(file.path(modelOutputDirectory,"ModelErrorRate.pdf"),width=6,height=6)
    thresholds=regularizationThresholds
    errorRates=thresholdCVRates[,"cvm"]
    ylim=c(0,1)
    if (family=="continuous"){
      ylim=c(min(thresholdCVRates[,"cvm"]-thresholdCVRates[,"cvsd"])*.9,max(thresholdCVRates[,"cvm"]+thresholdCVRates[,"cvsd"])*1.1)
    }
    ylab="Model Cross Validation Error Rate"
    if (modelType=="glmnet"){
      thresholds = log(thresholds)
      xlab="log(Regularization Threshold)"
      nonzeroCounts = finalModel$df
      if (family=="survival"){
        ylim=range(errorRates,na.rm=T)
        ylab="Model Cross Validation Partial Likelihood"
      }
    } else if (modelType=="pamr"){
      xlab="Regularization Threshold"
      nonzeroCounts = finalModel$nonzero
    }
    plot(errorRates,type='o',pch=20,col="red",main="Number of model features\n",axes=F,xlab=xlab,ylim=ylim,ylab=ylab)
    #Plot SEM
    for (i in 1:length(thresholds)){
      lines(c(i,i),c(errorRates[i]+thresholdCVRates[,"cvsd"][i],errorRates[i]-thresholdCVRates[,"cvsd"][i]),col="red",lty=3)
    }
    grid()
    axis(1,at=1:length(errorRates),labels=sapply(thresholds,.formatDecimal))
    if (family=="survival"){
      axis(2)
    } else {
      axis(2,at=c(0,.25,.5,.75,1),labels=c(0,25,50,75,100))  
    }
    
    axis(3,labels=nonzeroCounts,at=1:length(errorRates))
    legendLabels = c("Cross Validation Error Rate")
    legendColors = c("red")
    legendPchs = c(20)
    legendLty = c(1)
    legendPtCex=c(1)
    if (!is.null(thresholdCVRates$fdr)){
      lines(thresholdCVRates$fdr,type='o',pch=1,cex=1.5,col="blue")
      legendLabels = c(legendLabels,"Feature False Discovery Rate")
      legendColors = c(legendColors,"blue")
      legendPchs = c(legendPchs,1)
      legendLty = c(legendLty,1)
      legendPtCex=c(legendPtCex,1)
    }
    
    cv.min = cvMinima$cv.min.index
    cv.1se = cvMinima$cv.1se.index
    if (!is.null(cv.min)){
      points(c(cv.min,cv.min),y=c(errorRates[cv.min],errorRates[cv.min]),col="green",pch=20,cex=2)  
    }
    if (!is.null(cv.1se)){
      points(c(cv.1se,cv.1se),y=c(errorRates[cv.1se],errorRates[cv.1se]),col="orange",pch=9,cex=2)  
    }
        
    legendLabels = c(legendLabels,"cv.min","cv.1se")
    legendColors = c(legendColors,"green","orange")
    legendPchs = c(legendPchs,20,9)
    legendLty = c(legendLty,0,0)
    legendPtCex=c(legendPtCex,2,1.5)  
       
    if ( "cv.fdr.constrained" %in% names(cvMinima)) {
      
      cv.fdr.constrained = cvMinima$cv.fdr.constrained.index
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="yellow",pch=17,cex=1.5)
      points(c(cv.fdr.constrained,cv.fdr.constrained),y=c(errorRates[cv.fdr.constrained],errorRates[cv.fdr.constrained]),col="black",pch=2,cex=1.5)
      legendLabels = c(legendLabels,"cv.fdr.constrained")
      legendColors = c(legendColors,"yellow")
      legendPchs = c(legendPchs,17)
      legendLty = c(legendLty,0)
      legendPtCex=c(legendPtCex,1.5)  
      
    }
    legend(x="topleft",legendLabels,col=legendColors,pch=legendPchs,lty=legendLty,pt.cex=legendPtCex,cex=.8,bg="white")
    dev.off()   
}

citrus.plotModelDifferentialFeatures.classification = function(differentialFeatures,features,modelOutputDirectory,labels,...){
  for (cvPoint in names(differentialFeatures)){
    nonzeroFeatureNames = differentialFeatures[[cvPoint]][["features"]]

    # Write features to file for easy parsing
    write.table(features[,nonzeroFeatureNames],file=file.path(modelOutputDirectory,paste("features_",cvPoint,".csv",sep="")),quote=F,sep=",")

    melted = melt(data.frame(features[,nonzeroFeatureNames,drop=F],labels=labels,check.names=F),id.vars="labels")
    
    
    pdf(file.path(modelOutputDirectory,paste("features-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=4,height=length(nonzeroFeatureNames)*1.5)
    p <- ggplot(melted, aes(x=factor(labels), y=value)) 
    p = p + facet_wrap(~variable,ncol=1) + geom_boxplot(outlier.colour=rgb(0,0,0,0),colour=rgb(0,0,0,.3)) + geom_point(aes(color=factor(labels)),alpha=I(0.25),shape=19,size=I(2)) + coord_flip() +  theme_bw() + ylab("") + xlab("") + theme(legend.position = "none")
    if (any(grepl(pattern="abundance",nonzeroFeatureNames))){
      p  = p+scale_y_log10() + ylab("Log10 scale")
    }
    print(p)
    dev.off()
  }
  
}

citrus.plotModelDifferentialFeatures.continuous = function(differentialFeatures,features,modelOutputDirectory,labels,...){
  for (cvPoint in names(differentialFeatures)){
    nonzeroFeatureNames = differentialFeatures[[cvPoint]][["features"]]
    
    # Write features to file for easy parsing
    write.table(features[,nonzeroFeatureNames],file=file.path(modelOutputDirectory,paste("features_",cvPoint,".csv",sep="")),quote=F,sep=",")
    
    melted = melt(data.frame(features[,nonzeroFeatureNames,drop=F],labels=labels,check.names=F),id.vars="labels")
    
    pdf(file.path(modelOutputDirectory,paste("features-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep="")),width=4,height=length(nonzeroFeatureNames)*1.5)
    p <- ggplot(melted, aes(x=value, y=labels)) 
    p = p + facet_wrap(~variable,ncol=1) + geom_point(size=I(2)) + theme_bw() + ylab("") + xlab("") + theme(legend.position = "none")
    if (any(grepl(pattern="abundance",nonzeroFeatureNames))){
      p  = p+scale_x_log10() + xlab("Log10 scale")
    }
    print(p)
    dev.off()
  }
  
}

citrus.overlapDensityPlot = function(clusterDataList,backgroundData){
  combined = data.frame(check.names=F,check.rows=F)
  for (clusterName in names(clusterDataList)){
    combined=rbind(combined,cbind(melt(clusterDataList[[clusterName]],varnames=c("row","marker")),clusterId=clusterName,src="Cluster"))
  }
  p = ggplot(data=combined,aes(x=value, y = ..scaled..,fill=src)) + geom_density() + facet_grid(clusterId~marker,scales="free")+geom_density(data=cbind(melt(backgroundData,varnames=c("row","marker")),src="Background"))+theme_bw()+scale_fill_manual(values = c("Background" = rgb(.3,.3,1,.2), "Cluster" = rgb(1,.3,.3,.5)))
  p = p+theme(legend.position="bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_blank())
  p = p+labs(fill="Distribution:")
  print(p)
}

citrus.plotModelClusters = function(differentialFeatures,modelOutputDirectory,clusterAssignments,citrus.combinedFCSSet,clusteringColumns,...){
  for (cvPoint in names(differentialFeatures)){
    clusterIds = as.numeric(differentialFeatures[[cvPoint]][["clusters"]])
    outputFile = file.path(modelOutputDirectory,paste("clusters-",sub(pattern="\\.",replacement="_",x=cvPoint),".pdf",sep=""))
    citrus.plotClusters(clusterIds,clusterAssignments=clusterAssignments,citrus.combinedFCSSet,clusteringColumns,outputFile=outputFile,...)
  }
}

#' Plot cluster histograms
#' 
#' Plot expression of markers in cluster cells relative to all cells
#' 
#' @param clusterIds Vector of cluster IDs to plot
#' @param clusterAssignments List containing indicies of cells assigned to each cluster.
#' @param citrus.combinedFCSSet Combined FCS data that was clustered.
#' @param clusteringColumns Columns for which to plot distributions
#' @param conditions Vector of conditions clustering was performed on.
#' @param outputFile If not \code{NULL}, plot is written to \code{outputFile}.
#' @param ... Other parameters (ignored).
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Plot clusters
#' citrus.plotClusters(clusterIds=c(19998,19997),clusterAssignments=citrus.clustering$clusterMembership,citrus.combinedFCSSet,clusteringColumns)
citrus.plotClusters = function(clusterIds,clusterAssignments,citrus.combinedFCSSet,clusteringColumns,conditions=NULL,outputFile=NULL,...){
  
  data = citrus.combinedFCSSet$data
  
  if (!is.null(outputFile)){
    pdf(file=outputFile,width=(2.2*length(clusteringColumns)+2),height=(2*length(clusterIds)))  
  }
  clusterDataList = list();
  for (clusterId in sort(clusterIds)){
    if (length(clusterAssignments[[clusterId]])>2500){
      clusterDataList[[as.character(clusterId)]]=data[clusterAssignments[[clusterId]],clusteringColumns][sample(1:length(clusterAssignments[[clusterId]]),2500),]
    } else {
      clusterDataList[[as.character(clusterId)]]=data[clusterAssignments[[clusterId]],clusteringColumns]
    }
    
    colnames(clusterDataList[[as.character(clusterId)]])=.getDisplayNames(citrus.combinedFCSSet,clusteringColumns)
    
  }
  if (nrow(data)>2500){
    bgData = data[sample(1:nrow(data),2500),clusteringColumns]
  } else {
    bgData = data[,clusteringColumns]
  }
  colnames(bgData) = colnames(clusterDataList[[1]])
  citrus.overlapDensityPlot(clusterDataList=clusterDataList,backgroundData=bgData)
  if (!is.null(outputFile)){
    dev.off()  
  }
}

########################################
# Hierarchical Clustering plots
########################################

#' Create a graph object from clustering hierarchy
#' 
#' Create a graph object from clustering hierarchy that may be used by other plotting functions.
#' 
#' @param citrus.foldFeatureSet A \code{citrus.foldFeatureSet} object. Clusters for which features are calculated are included in the graph.
#' @param citrus.foldClustering A \code{citrus.foldClustering} object. Used to determine relationships between clusters.
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @return A \code{graph} object and other properties for plotting. 
#' \item{graph}{An \code{\link{igraph}} \code{graph} object.}
#' \item{layout}{An \code{\link{igraph}} \code{layout} for the graph.}
#' \item{plotSize}{Size of PDF for plotting (inches) calculated based on the number of clusters to be plotted.}
#' 
#' @examples
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Create graph for plotting
#' hierarchyGraphStuff = citrus.createHierarchyGraph(citrus.clustering,selectedClusters=largeEnoughClusters)
citrus.createHierarchyGraph = function(citrus.clustering,selectedClusters){
  
  if (length(selectedClusters)>750){
    minVertexSize=0
    plotSize=35
  } else if (length(selectedClusters)>250){
    minVertexSize=4
    plotSize=20
  } else if (length(selectedClusters)>100){
    minVertexSize=6
    plotSize=15
  } else {
    minVertexSize=8
    plotSize=10
  }
  
  
  mergeOrder=citrus.clustering$clustering$merge
  clusterAssignments=citrus.clustering$clusterMembership
  
  cmat = matrix(0,ncol=length(selectedClusters),nrow=length(selectedClusters),dimnames=list(selectedClusters,selectedClusters))
  for (cluster in selectedClusters){
    children = mergeOrder[cluster,]  
    for (child in children){
      if (as.character(child) %in% rownames(cmat)){
        cmat[as.character(cluster),as.character(child)] = 1
      }
    }
  }
  g = graph.adjacency(adjmatrix=cmat,mode="directed",add.colnames='label')
  clusterSizes = do.call("rbind",lapply(clusterAssignments,length))
  sizes = sapply(log(clusterSizes[selectedClusters]),findInterval,vec=hist(log(clusterSizes[selectedClusters]),breaks=10,plot=F)$breaks)
  sizes = sizes+minVertexSize
  g = set.vertex.attribute(g,"size",value=sizes)
  l = layout.reingold.tilford(g,root=length(V(g)),circular=T)
  result = list(graph=g,layout=l,plotSize=plotSize)  
  return(result)
}

#' Plot clustering hierarchy
#' 
#' Plots clustering hierarchy in graph form
#' 
#' @param outputFile Full path to output file (should have '.pdf' extension)
#' @param clusterColors Numeric or Character matrix of values to colors clusters by. See \code{details}.
#' @param graph Graph object to be plotted
#' @param layout Layout for graph
#' @param theme General color theme for plot. Options are \code{'black'} and \code{'white'}.
#' @param plotSize Size of square pdf (inches).
#' @param singlePDF Plot graphs for all variables in \code{clusterColors} in a single PDF?
#' @param ncol Number of columns if plotting all graphs in single PDF.
#' @param scale Scale up the size of the single PDF plot. 
#' @param plotClusterIDs Plot cluster IDs on vertices? 
#' 
#' @author Robert Bruggner
#' @export
#' 
#' @seealso \code{\link{citrus.createHierarchyGraph}}
#' 
#' @details The \code{clusterCols} argument enables multiple plots of the clustering hierarchy to be made, each colored 
#' by a different variable. \code{clusterCols} should be a numeric matrix with each cluster being plotted represented 
#' in a different row and each variable to be plotted represented in a different column. Row and column names should be 
#' cluster IDs and variable names respectively. If \code{clusterCols} is numeric, a color scale is generated across the
#' range of matrix values. Alternatively \code{clusterCols} can be a matrix of color names that are directly used to color
#' vertices.  
#' 
#' @examples
#' ############
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Create graph for plotting
#' hierarchyGraph = citrus.createHierarchyGraph(citrus.clustering,selectedClusters=largeEnoughClusters)
#' 
#' # Create matrix of variables to plot (in this case, cluster medians)
#' clusterMedians = t(sapply(largeEnoughClusters,citrus:::.getClusterMedians,clusterAssignments=citrus.clustering$clusterMembership,data=citrus.combinedFCSSet$data,clusterCols=clusteringColumns))
#' rownames(clusterMedians) = largeEnoughClusters
#' colnames(clusterMedians) = clusteringColumns
#' 
#' # Plot Clustering Hierarchy - Uncomment and Specify an output file
#' # citrus.plotClusteringHierarchy(outputFile="/path/to/output.pdf",clusterColors=clusterMedians,graph=hierarchyGraph$graph,layout=hierarchyGraph$layout,plotSize=hierarchyGraph$plotSize)
citrus.plotClusteringHierarchy = function(outputFile,clusterColors,graph,layout,theme="black",plotSize=15,singlePDF=F,ncol=3,scale=1,plotClusterIDs=T){
  if (theme=="black"){
    bg="black"
    stroke="white"
    strokea=rgb(1,1,1,.5)
  } else if (theme=="white"){
    bg="white"
    stroke="black"
    strokea=rgb(0,0,0,.5)
  } else {
    stop("Unrecognized theme option. Choices are 'white' or 'black'")
  }
  if (plotClusterIDs){
    vc="white"
  } else {
    vc=rgb(0,0,0,0)
  }
  if (singlePDF){
    expand=3*scale
    nrow = ceiling(ncol(clusterColors)/ncol)
    pdf(file=outputFile,width=ncol*expand,height=nrow*expand,bg=bg)
    par(mfrow=c(nrow,ncol),oma=rep(0.8,4),mar=c(.1,0,.8,2.8))
  } else {
    pdf(file=outputFile,width=plotSize,height=plotSize,bg=bg)  
    
  }
  
  for (target in 1:ncol(clusterColors)){
    if (is.numeric(clusterColors)){
      ct = seq(from=(min(clusterColors[,target])-0.01),to=(max(clusterColors[,target])+0.01),length.out=20)
      cols = .graphColorPalette(20)[sapply(clusterColors[,target],findInterval,vec=ct)]  
    } else {
      cols = clusterColors[,target]
    }
    par(col.main=stroke)  
    plot.igraph(graph,layout=layout,vertex.color=cols,main=colnames(clusterColors)[target],edge.color=stroke,vertex.label.color=vc,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica")
        
    # Legend
    if(is.numeric(clusterColors)){
      legend_image <- as.raster(matrix(rev(.graphColorPalette(20)), ncol=1))
      rasterImage(legend_image, 1.1, -.5, 1.15,.5)
      text(x=1.15, y = seq(-.5,.5,l=5), labels = .decimalFormat(ct[c(1,floor((length(ct)/4)*1:4))]) ,pos=4,col=stroke)  
    }
  }
  dev.off()  
}

#' Plot clustering hierarchy with clusters highlighted
#' 
#' Plot clustering heirarchy with a subset of clusters (clusters of interest) highlighted by color and/or encircled. 
#' 
#' @param outputFile Full path to output file (should end in '.pdf'). 
#' @param featureClusterMatrix Matrix of clusters to encircle. See details.
#' @param graph Graph object to be plotted
#' @param layout Layout for graph
#' @param theme General color theme for plot. Options are \code{'black'} and \code{'white'}.
#' @param plotSize Size of square pdf (inches).
#' @param plotClusterIDs Plot cluster IDs on vertices? 
#' @param featureClusterColors Named vector of colors for each vertex. See details.
#' @param encircle Should related highlighted clusters be encircled? 
#' 
#' @details The \code{featureClusterMatrix} argument should be a two column matrix. The first column should be names 'cluster' and rows should contain
#' a cluster id to be highlighted. The second column should be named 'feature' and should contain a string or property describing 
#' the general property of interest for this cluster. Entries having the same 'feature' value are plotted in the graph. One plot 
#' is created for each unique value of the 'feature' column. 
#' 
#' The \code{featureClusterColors} argument can be used to supply custom colors for graph vertices. The vector should be a named
#' vector with each entry being a color value and the name of the entry should be a vertex id (cluster id) that will be colored. 
#' 
#' @author Robert Bruggner
#' @export
#'
#' @examples 
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # Cluster data
#' citrus.clustering = citrus.cluster(citrus.combinedFCSSet,clusteringColumns)
#' 
#' # Large enough clusters
#' largeEnoughClusters = citrus.selectClusters(citrus.clustering)
#' 
#' # Create graph for plotting
#' hierarchyGraph = citrus.createHierarchyGraph(citrus.clustering,selectedClusters=largeEnoughClusters)
#' 
#' # Features to highlight 
#' featureClusterMatrix = data.frame(cluster=c(19992,19978,19981,19987,19983,19973),feature=rep(c("Property 1","Property 2"),each=3))
#' 
#' # Plot features in clustering hierarchy 
#' # citrus.plotHierarchicalClusterFeatureGroups(outputFile="/path/to/outputFile.pdf",featureClusterMatrix,graph=hierarchyGraph$graph,layout=hierarchyGraph$layout,plotSize=hierarchyGraph$plotSize)
citrus.plotHierarchicalClusterFeatureGroups = function(outputFile,featureClusterMatrix,graph,layout,theme="black",plotSize=15,plotClusterIDs=T,featureClusterColors=NULL,encircle=T){
  
  if (!is.null(featureClusterColors)&&(is.null(names(featureClusterColors)))){
    stop("featureClusterColors argument must be vector with elements having names of vertices to be colored.")
  }
  
  if (theme=="black"){
    bg="black"
    stroke="white"
    strokea=rgb(1,1,1,.5)
  } else if (theme=="white"){
    bg="white"
    stroke="black"
    strokea=rgb(0,0,0,.5)
  } else {
    stop("Unrecognized theme option. Choices are 'white' or 'black'")
  }
  pdf(file=outputFile,width=plotSize,height=plotSize,bg=bg)
  for (feature in unique(featureClusterMatrix[,"feature"])){
    fGroup = list();
    featureClusters = featureClusterMatrix[featureClusterMatrix[,"feature"]==feature,"cluster"]
    featureElements = match(featureClusters,as.numeric(get.vertex.attribute(graph,"label",V(graph))))
    subgraph = induced.subgraph(graph,featureElements)
    groupAssignments = clusters(subgraph)$membership
    for (groupId in unique(groupAssignments)){
      fGroup[[groupId]]=match(get.vertex.attribute(subgraph,"label")[groupAssignments==groupId],get.vertex.attribute(graph,"label"))
    }
    par(col.main=stroke)

    vertexFont=rep(1,length(V(graph)))
    vertexFont[get.vertex.attribute(graph,"label")%in%featureClusters]=2
    if (is.null(featureClusterColors)){
      vertexColor=rep(rgb(0,0,.5,.5),length(V(graph)))
      vertexColor[get.vertex.attribute(graph,"label")%in%featureClusters]=rgb(0.5,0,0,.7)
    } else {
      vertexColor=rep(rgb(0,0,.5,.5),length(V(graph)))
      cp = rgb(1,0,0,seq(0,1,by=.05))
      ct = seq(from=(min(featureClusterColors)-0.01),to=(max(featureClusterColors)+0.01),length.out=20)
      vertexColor[ match(names(featureClusterColors),get.vertex.attribute(graph,"label")) ] = cp[sapply(featureClusterColors,findInterval,vec=ct)]
    }
    
    vertexLabelColor="white"
    if (!plotClusterIDs){
      vertexLabelColor = rgb(0,0,0,0)
    }
    if (encircle){
      plot.igraph(graph,layout=layout,mark.groups=fGroup,mark.expand=5,main=feature,edge.color=stroke,vertex.label.color=vertexLabelColor,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=vertexColor,mark.border=strokea,mark.col=.graphColorPalette(length(fGroup),alpha=.2),vertex.label.font=vertexFont)      
    } else {
      plot.igraph(graph,layout=layout,main=feature,edge.color=stroke,vertex.label.color=vertexLabelColor,edge.arrow.size=.2,vertex.frame.color=strokea,vertex.label.cex=.7,vertex.label.family="Helvetica",vertex.color=vertexColor,vertex.label.font=vertexFont)     
    }
    if (!is.null(featureClusterColors)){
      legend_image <- as.raster(matrix(rev(cp), ncol=1))
      rasterImage(legend_image, 1.1, -.5, 1.15,.5)
      text(x=1.15, y = seq(-.5,.5,l=5), labels = .decimalFormat(ct[c(1,floor((length(ct)/4)*1:4))]) ,pos=4,col=stroke)  
    }
          
  }
  dev.off()
}

#' Plot results of a Citrus regression analysis
#' 
#' Makes many plots showing results of a Citrus analysis
#' 
#' @name plot.citrus.regressionResult
#' @export 
#' 
#' @param citrus.regressionResult A \code{citrus.regressionResult} object.
#' @param outputDirectory Full path to output directory for plots.
#' @param citrus.foldClustering A \code{citrus.foldClustering} object.
#' @param citrus.foldFeatureSet A \code{citrus.foldFeatureSet} object. 
#' @param citrus.combinedFCSSet A \code{citrus.combinedFCSSet} object.
#' @param plotTypes Vector of plots types to make. Valid options are \code{errorRate} (Cross-validated error rates for predictive models),
#' \code{stratifyingFeatures} (plots of non-zero model features),\code{stratifyingClusters} (plots of clustering marker distributions in stratifying clusters), and
#' \code{clusterGraph} (Plots of clustering hierarchy graph).
#' @param hierarchyGraph A hierarchy graph configuration created by \code{\link{citrus.createHierarchyGraph}}. If \code{NULL}, automatically generated.
#' 
#' @author Robert Bruggner
#' 
#' @examples 
#' # Where the data lives
#' dataDirectory = file.path(system.file(package = "citrus"),"extdata","example1")
#' 
#' # Create list of files to be analyzed
#' fileList = data.frame("unstim"=list.files(dataDirectory,pattern=".fcs"))
#' 
#' # Read the data 
#' citrus.combinedFCSSet = citrus.readFCSSet(dataDirectory,fileList)
#' 
#' # List of columns to be used for clustering
#' clusteringColumns = c("Red","Blue")
#' 
#' # List disease group of each sample
#' labels = factor(rep(c("Healthy","Diseased"),each=10))
#' 
#' # Cluster data
#' citrus.foldClustering = citrus.clusterAndMapFolds(citrus.combinedFCSSet,clusteringColumns,nFolds=1)
#' 
#' # Build abundance features
#' citrus.foldFeatureSet = citrus.calculateFoldFeatureSet(citrus.foldClustering,citrus.combinedFCSSet)
#' 
#' # Endpoint regress
#' citrus.regressionResult = citrus.endpointRegress(modelType="pamr",citrus.foldFeatureSet,labels,family="classification")
#' 
#' # Plot results
#' # plot(citrus.regressionResult,outputDirectory,"/path/to/output/directory/",citrus.foldClustering,citrus.foldFeatureSet,citrus.combinedFCSSet)
plot.citrus.regressionResult = function(citrus.regressionResult,outputDirectory,citrus.foldClustering,citrus.foldFeatureSet,citrus.combinedFCSSet,plotTypes=c("errorRate","stratifyingFeatures","stratifyingClusters","clusterGraph"),hierarchyGraph=NULL,...){
addtlArgs = list(...)
  
  theme="black"
  if ("theme" %in% names(addtlArgs)){
    theme = addtlArgs[["theme"]]
  }
    
  modelType=citrus.regressionResult$modelType
  
  cat(paste("Plotting results for",modelType,"\n"))
  
  # Make model output directoy
  modelOutputDirectory = file.path(outputDirectory,paste0(modelType,"_results"))
  dir.create(modelOutputDirectory,showWarnings=F,recursive=T)
  
  if ("errorRate" %in% plotTypes){
    cat("Plotting Error Rate\n")
    citrus.plotTypeErrorRate(modelType=modelType,modelOutputDirectory=modelOutputDirectory,regularizationThresholds=citrus.regressionResult$regularizationThresholds,thresholdCVRates=citrus.regressionResult$thresholdCVRates,finalModel=citrus.regressionResult$finalModel$model,cvMinima=citrus.regressionResult$cvMinima,family=citrus.regressionResult$family)
  }
  
  if ("stratifyingFeatures" %in% plotTypes){
    cat("Plotting Stratifying Features\n")
    do.call(paste("citrus.plotModelDifferentialFeatures",citrus.regressionResult$family,sep="."),args=list(differentialFeatures=citrus.regressionResult$differentialFeatures,features=citrus.foldFeatureSet$allFeatures,modelOutputDirectory=modelOutputDirectory,labels=citrus.regressionResult$labels))    
  }
  
  if ("stratifyingClusters" %in% plotTypes){
    cat("Plotting Stratifying Clusters\n")
    citrus.plotModelClusters(differentialFeatures=citrus.regressionResult$differentialFeatures,modelOutputDirectory=modelOutputDirectory,clusterAssignments=citrus.foldClustering$allClustering$clusterMembership,citrus.combinedFCSSet=citrus.combinedFCSSet,clusteringColumns=citrus.foldClustering$allClustering$clusteringColumns,...)
  }
  
  
  if ("clusterGraph" %in% plotTypes){
    cat("Plotting Clustering Hierarchy")
    
    # Configure hierarchy graph if not provided
    if (is.null(hierarchyGraph)){
      hierarchyGraph = citrus.createHierarchyGraph(citrus.clustering=citrus.foldClustering$allClustering,selectedClusters=citrus.foldFeatureSet$allLargeEnoughClusters)  
    }
    
    
    # Plot median of clusters
    clusterMedians = t(sapply(citrus.foldFeatureSet$allLargeEnoughClusters,.getClusterMedians,clusterAssignments=citrus.foldClustering$allClustering$clusterMembership,data=citrus.combinedFCSSet$data,clusterCols=citrus.foldClustering$allClustering$clusteringColumns))
    rownames(clusterMedians) = citrus.foldFeatureSet$allLargeEnoughClusters
    colnames(clusterMedians) = .getDisplayNames(citrus.combinedFCSSet,clusteringColumns)
    
    citrus.plotClusteringHierarchy(outputFile=file.path(outputDirectory,"markerPlots.pdf"),clusterColors=clusterMedians,graph=hierarchyGraph$graph,layout=hierarchyGraph$layout,plotSize=hierarchyGraph$plotSize,theme=theme)
    citrus.plotClusteringHierarchy(outputFile=file.path(outputDirectory,"markerPlotsAll.pdf"),clusterColors=clusterMedians,graph=hierarchyGraph$graph,layout=hierarchyGraph$layout,plotSize=hierarchyGraph$plotSize,theme=theme,singlePDF=T,plotClusterIDs=F)
      
    for (cvPoint in names(citrus.regressionResult$differentialFeatures)){
      featureClusterMatrix = .getClusterFeatureMatrix(citrus.regressionResult$differentialFeatures[[cvPoint]][["features"]])
      citrus.plotHierarchicalClusterFeatureGroups(outputFile=file.path(modelOutputDirectory,paste("featurePlots_",cvPoint,".pdf",sep="")),featureClusterMatrix=featureClusterMatrix,graph=hierarchyGraph$graph,layout=hierarchyGraph$layout,plotSize=hierarchyGraph$plotSize,theme=theme)
    }
  }
  
}

#' Plot a citrus.full.result
#' @name plot.citrus.full.result
#' @param citrus.full.result A citrus.full.result object
#' @param outputDirectory Full path to directory in which to place plot output. 
#' 
#' @export
#' 
#' @author Robert Bruggner
#' 
#' @details See \code{\link{citrus.full}} for examples.
plot.citrus.full.result = function(citrus.full.result,outputDirectory){
  
  if (!file.exists(outputDirectory)){
    stop(paste0("Output directory '",outputDirectory,"' does not exist."))
  }
  # Should 
  for (conditionName in names(results$conditions)){
    cat(paste0("\nPlotting Results for ",conditionName,"\n"))
    conditionOutputDir = file.path(outputDirectory,conditionName)
    dir.create(conditionOutputDir,showWarnings=F)
    mclapply(citrus.full.result$conditionRegressionResults[[conditionName]],plot,outputDirectory=conditionOutputDir,citrus.foldClustering=citrus.full.result$citrus.foldClustering,citrus.foldFeatureSet=citrus.full.result$conditionFoldFeatures[[conditionName]],citrus.combinedFCSSet=citrus.full.result$citrus.combinedFCSSet,family=citrus.full.result$family,labels=citrus.full.result$labels,conditions=citrus.full.result$conditions[[conditionName]])
    cat("\n")
  }
}



