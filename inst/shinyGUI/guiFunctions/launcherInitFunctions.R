
getClusterCols = function(fileName,dataDir){
  colnames(read.FCS(file.path(dataDir,fileName),which.lines=1))
}
