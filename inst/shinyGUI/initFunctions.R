
getClusterCols = function(fileName,dataDir){
  colnames(read.FCS(paste(dataDir,fileName,sep=""),which.lines=1))
}
